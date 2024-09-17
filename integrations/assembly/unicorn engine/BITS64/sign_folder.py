# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is how to come close to putting a linux binary inside unicorn engine I got.

With some polish and elbow grease, it will run perfectly fine, as long as all the "right pieces" are linked to the right place.

"""


from unicorn import *
from unicorn.x86_const import *
from capstone import *
from elftools.elf.elffile import ELFFile

# Define memory addresses and sizes
BASE_ADDRESS = 0x40000
STACK_ADDRESS = 0x7ff00000  # Stack address aligned to page boundary
STACK_SIZE = 20 * 1024 * 1024  # 20 MB stack size for broader coverage
MAP_REGION_SIZE = 0x1000000  # 16 MB mapped region size to cover more areas
INITIAL_STUB_ADDRESS = 0x200000  # Initial address for the stub function

def align_address(addr):
    """Align address to the nearest 4KB boundary."""
    return (addr + 0xFFF) & ~0xFFF

def map_large_region(uc, base_addr, size):
    """Map a large memory region to avoid multiple smaller mappings."""
    try:
        uc.mem_map(base_addr, size, UC_PROT_READ | UC_PROT_WRITE | UC_PROT_EXEC)
        print(f"Large memory region mapped at {hex(base_addr)} with size {size} bytes")
    except UcError as e:
        print(f"Failed to map large memory region: {e}")
        return False
    return True

def check_address_conflict(uc, address, size):
    """Check if the desired address range conflicts with existing mappings."""
    for start, end, _ in uc.mem_regions():
        if address < end and (address + size) > start:
            return True
    return False

def find_free_address(uc, start_address, size):
    """Find a free address for the stub function, avoiding conflicts."""
    address = start_address
    while check_address_conflict(uc, address, size):
        print(f"Address conflict detected at {hex(address)}. Trying a new address.")
        address += 0x1000  # Move to the next 4KB-aligned address
    return address

def load_elf(uc, elf_path):
    with open(elf_path, 'rb') as f:
        elf = ELFFile(f)

        # Pre-map a large enough region to hold all ELF segments and additional space
        if not map_large_region(uc, BASE_ADDRESS, MAP_REGION_SIZE):
            return None

        # Load the ELF segments into this large region
        for segment in elf.iter_segments():
            if segment['p_type'] == 'PT_LOAD':
                mem_size = segment['p_memsz']
                mem_addr = BASE_ADDRESS + segment['p_vaddr']  # Use BASE_ADDRESS as a starting point
                data = segment.data()

                # Align memory size to the nearest 4KB
                aligned_size = align_address(mem_size)

                # Determine permissions for the segment
                perms = 0
                if segment['p_flags'] & 1:  # Executable
                    perms |= UC_PROT_EXEC
                if segment['p_flags'] & 2:  # Writable
                    perms |= UC_PROT_WRITE
                if segment['p_flags'] & 4:  # Readable
                    perms |= UC_PROT_READ

                try:
                    uc.mem_write(mem_addr, data)
                    print(f"Loaded segment at {hex(mem_addr)} with size {aligned_size} bytes")
                except UcError as e:
                    print(f"Failed to load segment at {hex(mem_addr)}: {e}")

        # Set the entry point of the binary
        entry_point = BASE_ADDRESS + elf.header['e_entry']
        print(f"Entry point is at: {hex(entry_point)}")
        return entry_point

def hook_code(uc, address, size, user_data):
    # Hook to print each instruction
    md = Cs(CS_ARCH_X86, CS_MODE_64)
    code = uc.mem_read(address, size)
    for i in md.disasm(code, address):
        print(f"0x{i.address:x}:\t{i.mnemonic}\t{i.op_str}")

def hook_mem_unmapped(uc, access, address, size, value, user_data):
    # Hook to catch unmapped memory access attempts
    print(f"Unmapped memory access at {hex(address)}")

    # Attempt to map additional memory dynamically
    aligned_address = align_address(address & ~0xFFF)  # Align to page boundary
    try:
        uc.mem_map(aligned_address, 0x1000, UC_PROT_READ | UC_PROT_WRITE | UC_PROT_EXEC)
        print(f"Dynamically mapped memory at {hex(aligned_address)} to handle access")
    except UcError as e:
        print(f"Failed to dynamically map memory at {hex(aligned_address)}: {e}")

def initialize_got_and_plt(uc, elf):
    # Identify GOT and PLT sections
    got_address = None
    plt_address = None
    got_size = 0
    plt_size = 0

    for section in elf.iter_sections():
        if section.name == '.got.plt' or section.name == '.got':
            got_address = BASE_ADDRESS + section['sh_addr']
            got_size = section['sh_size']
        elif section.name == '.plt':
            plt_address = BASE_ADDRESS + section['sh_addr']
            plt_size = section['sh_size']

    if got_address and plt_address:
        print(f"Initializing GOT at {hex(got_address)} with size {got_size} bytes.")
        print(f"Initializing PLT at {hex(plt_address)} with size {plt_size} bytes.")

        # Find a conflict-free address for the stub function
        stub_address = find_free_address(uc, INITIAL_STUB_ADDRESS, 0x1000)

        # Initialize GOT entries with the address of the stub function
        for offset in range(0, got_size, 8):  # Assuming 64-bit GOT entries
            entry_address = got_address + offset
            uc.mem_write(entry_address, stub_address.to_bytes(8, byteorder='little'))
            print(f"Set GOT entry at {hex(entry_address)} to {hex(stub_address)}")

        # Stub function code to handle unresolved function calls safely
        stub_code = b'\x48\x31\xc0\x48\xff\xc0\xc3'  # x86_64: xor rax, rax; inc rax; ret
        try:
            uc.mem_map(stub_address, 0x1000, UC_PROT_READ | UC_PROT_WRITE | UC_PROT_EXEC)
            uc.mem_write(stub_address, stub_code)
            print(f"Stub function mapped at {hex(stub_address)} to handle unresolved calls.")
        except UcError as e:
            print(f"Failed to map stub function at {hex(stub_address)}: {e}")

def main():
    # Initialize Unicorn Engine for x86_64
    uc = Uc(UC_ARCH_X86, UC_MODE_64)

    # Map a stack for the binary
    try:
        uc.mem_map(STACK_ADDRESS, STACK_SIZE, UC_PROT_READ | UC_PROT_WRITE)
        uc.reg_write(UC_X86_REG_RSP, STACK_ADDRESS + STACK_SIZE - 8)  # Set stack pointer near the top of the stack
        print(f"Stack mapped at {hex(STACK_ADDRESS)} with size {STACK_SIZE} bytes")
    except UcError as e:
        print(f"Failed to map stack: {e}")
        return

    # Load the ELF file
    elf_path = "/bin/ls"  # Replace with your ELF file path
    with open(elf_path, 'rb') as f:
        elf = ELFFile(f)
        entry_point = load_elf(uc, elf_path)

        if entry_point is None:
            return  # Exit if loading the ELF file failed

        # Initialize GOT and PLT if the ELF is dynamically linked
        initialize_got_and_plt(uc, elf)

    # Set a code hook to print instructions being executed
    uc.hook_add(UC_HOOK_CODE, hook_code)

    # Set a hook for unmapped memory fetch to catch and handle dynamically
    uc.hook_add(UC_HOOK_MEM_FETCH_UNMAPPED, hook_mem_unmapped)

    # Start emulating from the entry point
    try:
        # Set the end address to prevent running past mapped regions
        end_address = BASE_ADDRESS + MAP_REGION_SIZE  # Ensure it stays within mapped region
        uc.emu_start(entry_point, end_address)
    except UcError as e:
        print(f"Unicorn Error during emulation: {e}")

if __name__ == "__main__":
    main()

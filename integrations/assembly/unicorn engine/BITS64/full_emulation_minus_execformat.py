# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This, with the right interception of stack, can be used for emulating the arch your app should use. In this case, it is meant to mimick x86-64 CPU arch of a stack setting up a basic 
 interrupt to the linux kernel.


"""


from unicorn import *
from unicorn.x86_const import *

# Typical high address for stack in a 64-bit Linux process
STACK_ADDR = 0x2000
STACK_SIZE = 4 * 1024  # 2MB stack size

# Initialize emulator in X86_64 mode
mu = Uc(UC_ARCH_X86, UC_MODE_64)

# Map memory for the stack
mu.mem_map(STACK_ADDR - STACK_SIZE, STACK_SIZE)

# Set the stack pointer
mu.reg_write(UC_X86_REG_RSP, STACK_ADDR)

# Example code to emulate
ADDRESS = 0x0000
CODE = (
    b"\x55"                      # push   %rbp
    b"\x48\x89\xe5"              # mov    %rsp, %rbp
    b"\xb8\x01\x00\x00\x00"      # mov    $0x1, %eax
    b"\xcd\x80"                  # int    0x80 (simulated syscall)
    b"\x5d"                      # pop    %rbp
)

# Map memory for the code
mu.mem_map(ADDRESS, 4 * 1024)
mu.mem_write(ADDRESS, CODE)

# Add hook to handle interrupts
def hook_intr(uc, intno, user_data):
    if intno == 0x80:
        print("Intercepted int 0x80 - Simulated syscall")
        #uc.reg_write(UC_X86_REG_EAX, 0)  # Set return value for syscall

mu.hook_add(UC_HOOK_INTR, hook_intr)


# Hook to handle unmapped memory accesses
def hook_unmapped(uc, access, address, size, value, user_data):
    if access == UC_MEM_WRITE:
        print(f">>> Unmapped memory write at 0x{address:x}, size = {size}, value = 0x{value:x}")
    else:  # UC_MEM_READ or UC_MEM_FETCH
        print(f">>> Unmapped memory read at 0x{address:x}, size = {size}")
        # Align the address to the nearest page boundary
        aligned_address = address & ~0xfff
        mem_size = 0x1000  # 4KB page size
        uc.mem_map(aligned_address, mem_size)
        print(f">>> Mapped 0x{mem_size:x} bytes at 0x{aligned_address:x}, {value}")
        sleep(10)
    return True  # Return True to indicate we handled this exception

mu.hook_add(UC_HOOK_MEM_UNMAPPED, hook_unmapped)

# Hook to print instructions
def hook_code(mu, address, size, user_data):
    print(f">>> Tracing instruction at 0x{address:x}, instruction size = {size}")
    # read the memory bytes that are being executed
    bytes = mu.mem_read(address, size)
    print(f">>> Instructions: {bytes.hex()}")

mu.hook_add(UC_HOOK_CODE, hook_code)


# Start emulation
try:
    mu.emu_start(ADDRESS, ADDRESS + len(CODE))
except UcError as e:
    print(f"ERROR: {e}")

# Print final register values
print(f"RAX = 0x{mu.reg_read(UC_X86_REG_RAX):x}")
print(f"RSP = 0x{mu.reg_read(UC_X86_REG_RSP):x}")

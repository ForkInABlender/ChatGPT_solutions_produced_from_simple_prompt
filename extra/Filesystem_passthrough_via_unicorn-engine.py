# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""

The following would be useful for filesystem emulation
 during passthrough. This also makes it easier to
  only use the minimum logic needed to run dockerd.
 This saves work as one can make use of a specific
  directory, and not have to worry about the runtime
 not saving where it last was ehen it crashed last.

"""

import os
import struct
from unicorn import *
from unicorn.x86_const import *

# Constants for system call numbers
SYS_open = 5
SYS_read = 3
SYS_write = 4
SYS_close = 6
SYS_lseek = 19

# Directory to emulate from
BASE_DIR = "/tmp/uc_emulator"

# Dictionary to map emulated file descriptors to real file descriptors
fd_map = {}

# Ensure the base directory exists
os.makedirs(BASE_DIR, exist_ok=True)

# Initialize Unicorn Engine
def init_unicorn():
    uc = Uc(UC_ARCH_X86, UC_MODE_32)
    # Map 2MB memory for the emulation
    ADDRESS = 0x1000000
    SIZE = 2 * 1024 * 1024
    uc.mem_map(ADDRESS, SIZE)
    return uc, ADDRESS

# Callback function for handling file-related system calls
def hook_intr(uc, intno, user_data):
    if intno == 0x80:
        eax = uc.reg_read(UC_X86_REG_EAX)
        ebx = uc.reg_read(UC_X86_REG_EBX)
        ecx = uc.reg_read(UC_X86_REG_ECX)
        edx = uc.reg_read(UC_X86_REG_EDX)

        if eax == SYS_open:
            filename = uc.mem_read(ebx, 256).split(b'\x00', 1)[0].decode('utf-8')
            flags = ecx
            real_path = os.path.join(BASE_DIR, filename)
            fd = os.open(real_path, flags)
            emu_fd = fd + 100  # Use an offset to avoid conflict with real FDs
            fd_map[emu_fd] = fd
            uc.reg_write(UC_X86_REG_EAX, emu_fd)
            print(f"open('{filename}', {flags}) = {emu_fd}")
        elif eax == SYS_read:
            fd = fd_map.get(ebx, -1)
            buffer_addr = ecx
            length = edx
            if fd != -1:
                data = os.read(fd, length)
                uc.mem_write(buffer_addr, data)
                uc.reg_write(UC_X86_REG_EAX, len(data))
                print(f"read({fd}, {length}) = {len(data)} bytes")
            else:
                uc.reg_write(UC_X86_REG_EAX, -1)
        elif eax == SYS_write:
            fd = fd_map.get(ebx, -1)
            buffer_addr = ecx
            length = edx
            if fd != -1:
                data = uc.mem_read(buffer_addr, length)
                written = os.write(fd, data)
                uc.reg_write(UC_X86_REG_EAX, written)
                print(f"write({fd}, {length}) = {written} bytes")
            else:
                uc.reg_write(UC_X86_REG_EAX, -1)
        elif eax == SYS_close:
            fd = fd_map.pop(ebx, -1)
            if fd != -1:
                os.close(fd)
                uc.reg_write(UC_X86_REG_EAX, 0)
                print(f"close({fd}) = 0")
            else:
                uc.reg_write(UC_X86_REG_EAX, -1)
        elif eax == SYS_lseek:
            fd = fd_map.get(ebx, -1)
            offset = ecx
            whence = edx
            if fd != -1:
                new_pos = os.lseek(fd, offset, whence)
                uc.reg_write(UC_X86_REG_EAX, new_pos)
                print(f"lseek({fd}, {offset}, {whence}) = {new_pos}")
            else:
                uc.reg_write(UC_X86_REG_EAX, -1)
        else:
            print(f"Unhandled system call: {eax}")

# Write machine code to simulate file operations
def write_machine_code(uc, address):
    # Machine code to open, read, write, and close a file
    X86_CODE32 = (
        # open("testfile", O_RDWR | O_CREAT, 0644)
        b"\xb8\x05\x00\x00\x00"  # mov eax, SYS_open
        + b"\xbb\x00\x20\x10\x00"  # mov ebx, filename address
        + b"\xb9\x42\x00\x00\x00"  # mov ecx, O_RDWR | O_CREAT
        + b"\xba\xb6\x01\x00\x00"  # mov edx, 0644
        + b"\xcd\x80"  # int 0x80
        
        # read(fd, buffer, 12)
        + b"\x89\xc3"  # mov ebx, eax (save fd)
        + b"\xb8\x03\x00\x00\x00"  # mov eax, SYS_read
        + b"\x89\xdf"  # mov edi, ebx (restore fd)
        + b"\xbb\x10\x30\x00\x00"  # mov ebx, buffer address
        + b"\xb9\x0c\x00\x00\x00"  # mov ecx, 12
        + b"\xcd\x80"  # int 0x80
        
        # write(fd, buffer, 12)
        + b"\xb8\x04\x00\x00\x00"  # mov eax, SYS_write
        + b"\x89\xdf"  # mov edi, ebx (restore fd)
        + b"\xbb\x10\x30\x00\x00"  # mov ebx, buffer address
        + b"\xb9\x0c\x00\x00\x00"  # mov ecx, 12
        + b"\xcd\x80"  # int 0x80
        
        # close(fd)
        + b"\xb8\x06\x00\x00\x00"  # mov eax, SYS_close
        + b"\x89\xdf"  # mov edi, ebx (restore fd)
        + b"\xcd\x80"  # int 0x80
    )
    # Write the machine code to memory
    uc.mem_write(address, X86_CODE32)
    # Write filename and buffer to memory
    filename = b"testfile\x00"
    buffer = b"Hello, World"
    uc.mem_write(address + 0x1020, filename)
    uc.mem_write(address + 0x1030, buffer)

# Main function to run the emulation
def main():
    # Initialize Unicorn Engine
    uc, address = init_unicorn()
    
    # Write machine code to simulate file operations
    write_machine_code(uc, address)
    
    # Set up the interrupt hook
    uc.hook_add(UC_HOOK_INTR, hook_intr)
    
    # Emulate the code
    try:
        uc.emu_start(address, address + 0x40)  # Adjust the length as needed
    except UcError as e:
        print(f"Unicorn error: {e}")
    
    print("Emulation done.")

# Run the main function
if __name__ == "__main__":
    main()

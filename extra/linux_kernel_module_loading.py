# Dylan Kenneth Eliot & GPT-4 ( Apha Edition )

"""

This is a way to load the module one would need for a application like dockerd.
 That's if you're trying to get it to run on android as python code without
  modding android ROM code.









"""

from unicorn import *
from unicorn.x86_const import *
import ctypes
import struct

# System call number for init_module
__NR_init_module = 128

# Callback function for handling system calls
def hook_intr(uc, intno, user_data):
    if intno == 0x80:
        eax = uc.reg_read(UC_X86_REG_EAX)
        ebx = uc.reg_read(UC_X86_REG_EBX)
        ecx = uc.reg_read(UC_X86_REG_ECX)
        edx = uc.reg_read(UC_X86_REG_EDX)

        if eax == __NR_init_module:
            print(f"System call: init_module(module_image=0x{ebx:x}, len={ecx}, param_values=0x{edx:x})")
            module_image = uc.mem_read(ebx, ecx)
            print(f"Module content (first 64 bytes): {module_image[:64].hex()}...")
            # Simulate module load success
            uc.reg_write(UC_X86_REG_EAX, 0)
        else:
            print(f"Unhandled system call: {eax}")

# Initialize emulator in X86-32bit mode
uc = Uc(UC_ARCH_X86, UC_MODE_32)

# Map 2MB memory for this emulation
ADDRESS = 0x1000000
SIZE = 2 * 1024 * 1024
uc.mem_map(ADDRESS, SIZE)

# Allocate space for the .ko file in the emulated memory
ko_address = ADDRESS + 0x1000
with open(ko_file_path, 'rb') as f:
    ko_content = f.read()
ko_size = len(ko_content)
uc.mem_write(ko_address, ko_content)

# Machine code to call init_module
X86_CODE32 = (
    b"\xb8\x80\x00\x00\x00"  # mov eax, 128 (__NR_init_module)
    + b"\xbb" + struct.pack("<I", ko_address)  # mov ebx, ko_address
    + b"\xb9" + struct.pack("<I", ko_size)     # mov ecx, ko_size
    + b"\xba\x00\x00\x00\x00"  # mov edx, 0 (param_values = NULL)
    + b"\xcd\x80"  # int 0x80
)

# Write the machine code to memory
uc.mem_write(ADDRESS, X86_CODE32)

# Set up the interrupt hook
uc.hook_add(UC_HOOK_INTR, hook_intr)

# Emulate the code
try:
    uc.emu_start(ADDRESS, ADDRESS + len(X86_CODE32))
except UcError as e:
    print(f"Unicorn error: {e}")

print("Emulation done.")

# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""


This is assembly and python code being used to make, turn on and off swap.

one edit later, and it should be enough to also delete the swap file upon termination.


"""

import subprocess
import os
import tempfile
from ctypes import *

PAGE_SIZE = 4096    

class AssemblerFunction(object):
    def __init__(self, code, ret_type, *arg_types):
        # Run NASM
        fd, source = tempfile.mkstemp(".asm", "assembly", os.getcwd())
        with os.fdopen(fd, 'w') as f:
            f.write(code)
        target = os.path.splitext(source)[0] + ".o"
        subprocess.check_call(["nasm", "-f", "elf64", source, "-o", target])
        os.unlink(source)
        binary = open(target, "rb").read()
        os.unlink(target)

        # Align our code on a page boundary
        self.code_buffer = create_string_buffer(PAGE_SIZE * 2 + len(binary))
        addr = (addressof(self.code_buffer) + PAGE_SIZE) & (~(PAGE_SIZE - 1))
        memmove(addr, binary, len(binary))

        # Change memory protection to executable
        self.mprotect = cdll.LoadLibrary("libc.so.6").mprotect
        if self.mprotect(addr, len(binary), 0x1 | 0x4):  # PROT_READ | PROT_EXEC
            raise OSError("Unable to change memory protection")

        self.func = CFUNCTYPE(ret_type, *arg_types)(addr)

    def __call__(self, *args):
        return self.func(*args)

    def __del__(self):
        # Revert memory protection
        if hasattr(self, "mprotect"):
            self.mprotect(self.addr, self.bin_len, 0x1 | 0x2)  # PROT_READ | PROT_WRITE

# Define the assembly code for each operation
swapon_code = """
section .data
swapfile db '/path/to/swapfile', 0

section .text
global _start
_start:
    mov eax, 167            ; swapon syscall number on Linux x86_64
    mov rdi, swapfile       ; argument: pointer to path
    xor esi, esi            ; flags, set to 0
    syscall                 ; make the syscall
    ret
"""

swapoff_code = """
section .data
swapfile db '/path/to/swapfile', 0

section .text
global _start
_start:
    mov eax, 168            ; swapoff syscall number on Linux x86_64
    mov rdi, swapfile       ; argument: pointer to path
    syscall                 ; make the syscall
    ret
"""

mkswap_code = """
section .data
swapfile db '/path/to/swapfile', 0
swap_signature db "SWAPSPACE2", 0
header_size equ 4096  ; assuming a header size of 4096 bytes for simplicity

section .bss
header resb header_size  ; reserve space for the header

section .text
global _start
_start:
    ; Open the swap file
    mov eax, 2            ; sys_open
    mov rdi, swapfile     ; filename
    mov rsi, 02h          ; flags (O_RDWR)
    xor rdx, rdx          ; mode (not relevant for opening existing files)
    syscall               ; open file
    mov ebx, eax          ; save file descriptor

    ; Prepare the header
    mov eax, [header]
    mov dword [eax], swap_signature ; Write the swap signature at start
    mov eax, [header + 10h]         ; Assume offset 10h for size
    mov dword [eax], 1024*1024*1024 ; Size of swap area (1GB for example)

    ; Write the header to the file
    mov eax, 1            ; sys_write
    mov edi, ebx          ; file descriptor
    mov rsi, header       ; buffer
    mov edx, header_size  ; number of bytes to write
    syscall

    ; Close the file
    mov eax, 3            ; sys_close
    mov edi, ebx          ; file descriptor
    syscall

    ; Exit
    mov eax, 60           ; sys_exit
    xor edi, edi          ; status 0
    syscall
"""

# Create instances of the AssemblerFunction for each operation
Swapon = AssemblerFunction(swapon_code, c_int)
Swapoff = AssemblerFunction(swapoff_code, c_int)
MkSwap = AssemblerFunction(mkswap_code, c_int)

# Example usage
if __name__ == "__main__":
    print("Setting up swap...")
    MkSwap()  # Initialize swap
    Swapon()
    input("press enter to turn off swap file.. :")
    Swapoff()
    print("swapfile off")

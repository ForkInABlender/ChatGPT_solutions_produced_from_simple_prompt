import subprocess
from ctypes import CFUNCTYPE, c_char_p, c_int, create_string_buffer

# Assembly code to be compiled
asm_code = """
section .data
    msg db 'Hello, World!', 0

section .text
    global _start

_start:
    ; write the message to stdout
    mov eax, 4 ; system call number (sys_write)
    mov ebx, 1 ; file descriptor (stdout)
    mov ecx, msg ; pointer to message
    mov edx, 13 ; length of message
    int 0x80 ; call kernel

    ; exit
    mov eax, 1 ; system call number (sys_exit)
    xor ebx, ebx ; exit status
    int 0x80 ; call kernel
"""

# link the assembly code with nasm
proc = subprocess.run(["nasm", "-f", "elf", "-o", "-"], input=asm_code, stdout=subprocess.PIPE)

# Load the binary into Python
from ctypes import CDLL

hello = CDLL(None, handle=proc.stdout)
hello._start.restype = c_int

# Access the msg data in the .data section
msg = c_char_p.in_dll(hello, "msg")
print(msg.value.decode())

# Create a Python function that wraps the _start function
@CFUNCTYPE(c_int)
def hello_function():
    return hello._start()

# Call the Python function
hello_function()

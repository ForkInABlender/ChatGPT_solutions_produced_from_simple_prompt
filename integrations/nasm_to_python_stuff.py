# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition) (code tune up)


"""
It is now updated to embed any assembly used directly into python. 


This can now directly bind to the assembly of any c/c++ binary that follows the elf binary format inside the python code with no disk bloatation.

On top of that, one likely could disassemble even docker and kubernetes tools and use their assembly directly inside python.

"""


import subprocess
import os
import tempfile
from ctypes import *

asm_code = """
section .data
    global msg
    msg db 'Hello, World!', 0xa

section .text
    global _start

_start:
    mov rax, 1
    mov rdi, 1
    mov rsi, msg
    mov rdx, 14
    syscall
   ret
"""

fd, asm_file = tempfile.mkstemp(suffix=".asm")
os.write(fd, asm_code.encode())
os.close(fd)
obj_file = asm_file.replace('.asm', '.o')
subprocess.run(["nasm", "-f", "elf64", asm_file, "-o", obj_file])
shared_lib = asm_file.replace('.asm', '.so')
subprocess.run(["ld", "-shared", "-o", shared_lib, obj_file])
hello = CDLL(shared_lib)
hello._start.argtypes = []
hello._start.restype = c_int
msg = (c_char*14).in_dll(hello, "msg")
os.remove(asm_file)
os.remove(obj_file)
os.remove(shared_lib)

hello._start()
print("normal")

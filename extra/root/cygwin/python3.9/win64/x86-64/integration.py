# Dylan Kenneth Eliot & GPT-4o-3mini0 ( alpha release on alpine using zipfile mounted as read only fusev3 )

"""
Enjoy. 

This now makes it easier to integrate other bits of interpreter code and decompilation to conform operating systems, kernels, threads, and processes, to make use of unicorn-engine
 where possible to minimize kernel exploitation vectors.


Thank you for your patience
  -- G-Null

"""


import os
import ctypes
import subprocess

# Write NASM assembly code to a temporary file
nasm_file = "init_python.asm"
so_file = "init_python.dll"

nasm_code = """
section .text
    global Py_InitializeEx
    extern Py_Initialize

Py_InitializeEx:
    push    rbp
    mov     rbp, rsp
    mov     edi, 1       ; Pass 1 to Py_Initialize (no signal handlers)
    call    Py_Initialize
    leave
    ret
"""

# Write NASM code to a file
with open(nasm_file, "w") as f:
    f.write(nasm_code)

# Assemble with NASM
subprocess.run(["nasm", "-f", "win64", nasm_file, "-o", "init_python.o"], check=True)

# Link with GCC
subprocess.run(["gcc", "-shared", "-o", so_file, "init_python.o", "-lpython3.9"], check=True)

# Load shared library
lib = ctypes.CDLL(so_file)
Py_InitializeEx = ctypes.CFUNCTYPE(None)(("Py_InitializeEx", lib))

# Call function
Py_InitializeEx()

print("✅ Python initialized successfully from NASM (Cygwin, Python 3.9)!")

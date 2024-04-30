# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition) (code tune up)

"""
This script, like the nasm version, uses assembly for the
 system it is running on.


Due note that it will take some time to manually define
 functionality the long and drawn out way. In the end,
  it should work enough to put dockerd the long way
 onto an android device.

No, there will not be support for Apple supported
 hardware or software, much like Microsoft as they're
  incapable of being open source software supporters.

"""

import subprocess
import os
import tempfile
from ctypes import *

# AArch64 Assembly code
asm_code = """
    .section .data
    .global msg
msg:
    .ascii "Hello, World!\\n"

    .section .text
    .global _start

_start:
    mov x0, #1            // file descriptor 1
    ldr x1, =msg          // pointer to msg
    mov x2, #14           // message length
    mov x8, #64           // syscall number (sys_write) for aarch64
    svc #0                // make syscall
    mov x0, #0            // status code
    mov x8, #93           // syscall number (sys_exit) for aarch64
    svc #0                // make syscall
"""

fd, asm_file = tempfile.mkstemp(suffix=".s")
os.write(fd, asm_code.encode())
os.close(fd)
obj_file = asm_file.replace('.s', '.o')
subprocess.run(["as", "-o", obj_file, asm_file])
shared_lib = asm_file.replace('.s', '.so')
subprocess.run(["ld", "-shared", "-o", shared_lib, obj_file])
hello = CDLL(shared_lib)
hello._start.argtypes = []
hello._start.restype = c_int
msg = (c_char * 14).in_dll(hello, "msg")
os.remove(asm_file)
os.remove(obj_file)
os.remove(shared_lib)

hello._start()
print("normal")

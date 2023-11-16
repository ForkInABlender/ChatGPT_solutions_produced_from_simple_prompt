# Dylan Kenneth Eliot & Google Bard AI

"""

What it does is compile golang in a separate process and then allow for binding to the functions directly from memory.

the  ``-buildmode=c-shared`` flag ensures that the compiled code is compatible with C/C++ code before it gets used in
 python.

"""


import subprocess
import ctypes

go_code = """
package main

import (
    "C"
    "fmt"
)

//export start
func start() {
    fmt.Println("Hello, World!")
}
"""
proc = subprocess.run(["go", "tool", "asm", "-o", "-"], input=go_code.encode("utf-8"), stdout=subprocess.PIPE)
hello = ctypes.CDLL(proc.stdout)
hello.start.restype = ctypes.c_void_p
hello.start()

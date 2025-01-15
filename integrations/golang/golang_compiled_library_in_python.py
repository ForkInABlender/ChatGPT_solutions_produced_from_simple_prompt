# Dylan Kenneth Eliot & Google Bard AI

"""

What it does is compile golang in a separate process and then allow for binding to the functions directly from memory.

the  ``-buildmode=c-shared`` flag ensures that the compiled code is compatible with C/C++ code before it gets used in
 python.


checked to work within userland.apk on android using jammy Ubuntu; compiles and runs without issue. -- 01/15/2025 @ 11:25 am



"""


import subprocess
import ctypes
import tempfile
import os

# Go code as a string
go_code = """
package main

import "C"

import (
	"fmt"
)

//export HelloWorld
func HelloWorld() {
	fmt.Println("Hello from Go!")
}

func main() {}
"""

temp_dir = tempfile.mkdtemp()
with tempfile.NamedTemporaryFile(dir=temp_dir, mode="w", delete=False, suffix=".go") as go_file:
    go_file.write(go_code)
    go_file_path = go_file.name
so_file_path = go_file_path.replace(".go", ".so")                                        # if in userland.apk:
subprocess.run(["go", "build", "-o", so_file_path, "-buildmode=c-shared", go_file_path]) # update to /usr/local/go/bin/go
os.remove(go_file_path)
hello = ctypes.CDLL(so_file_path)
os.remove(so_file_path)
hello.HelloWorld.argtypes = []
hello.HelloWorld.restype = ctypes.c_void_p
hello.HelloWorld()

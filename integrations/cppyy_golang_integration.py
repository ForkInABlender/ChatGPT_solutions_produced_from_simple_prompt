# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )


"""

This combines golang and c++ inside python code enough to allow separation of concerns. Even for python reannotated functions.
 This also allows for multiple parts to run with as little compiling as possible.



"""


import cppyy
import subprocess
import ctypes
import tempfile
import os

# Step 1: Go Code
go_code = """
package main

import "C"
import "fmt"

//export HelloWorld
func HelloWorld() {
    fmt.Println("Hello from Go!")
}

func main() {}
"""

# Create a temporary directory and Go file
temp_dir = tempfile.mkdtemp()
with tempfile.NamedTemporaryFile(dir=temp_dir, mode="w", delete=False, suffix=".go") as go_file:
    go_file.write(go_code)
    go_file_path = go_file.name

# Compile the Go code to a shared object file
so_file_path = go_file_path.replace(".go", ".so")
subprocess.run(["go", "build", "-o", so_file_path, "-buildmode=c-shared", go_file_path])

# Load the shared object file using ctypes
hello_lib = ctypes.CDLL(so_file_path)
hello_func = hello_lib.HelloWorld
hello_func.argtypes = []
hello_func.restype = None

# Step 2: Define the C++ Code
cpp_code = """
extern "C" {
    typedef void (*func_ptr)();

    void set_callback(func_ptr f);
    void call_callback();
}

func_ptr global_callback;

void set_callback(func_ptr f) {
    global_callback = f;
}

void call_callback() {
    if (global_callback) {
        global_callback();
    }
}
"""

# Include the C++ code using cppyy
cppyy.cppdef(cpp_code)

# Step 3: Convert the Go function to a C function pointer
CALLBACK_FUNC = ctypes.CFUNCTYPE(None)
c_callback = CALLBACK_FUNC(hello_func)

# Step 4: Set the callback in the C++ code and call it
cppyy.gbl.set_callback(c_callback)
cppyy.gbl.call_callback()

# Clean up temporary files and directories
os.remove(go_file_path)
os.remove(so_file_path)
os.rmdir(temp_dir)

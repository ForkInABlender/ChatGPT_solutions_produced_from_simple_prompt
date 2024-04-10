# Dylan Kenneth Eliot & Google Bard AI & GPT-4-plugins (Alpha Edition)

"""

Does similar to the c code in python script but for c++ code. Further extending the languages.
Now it can be used directly from python without needing something like cppyy.


"""


import subprocess
import ctypes
import tempfile
import os

cpp_code = """
#include <iostream>

    using namespace std;

    void start() {
        cout << "Hello, World!" << endl;
    }
"""

with tempfile.NamedTemporaryFile(suffix='.cpp', mode='w+', delete=False) as src_file:
    src_file.write(cpp_code)
    src_file_path = src_file.name

with tempfile.NamedTemporaryFile(suffix='.so', delete=False) as so_file:
    shared_object_path = so_file.name

subprocess.run(["g++", "-shared", "-fPIC", src_file_path, "-o", shared_object_path], check=True)

hello = ctypes.CDLL(shared_object_path)

os.system("objdump -T "+shared_object_path) #how to find the symbols to call
#hello.start() # <----- what is supposed to be available
hello._Z5startv() # <----- What the above ``os.system`` will be available until rewrapped with cython types in pure python mode..

os.remove(src_file_path)
os.remove(shared_object_path)

#Remove bloat to prevent excessive clutter. We don't want or need bloat to stay around.

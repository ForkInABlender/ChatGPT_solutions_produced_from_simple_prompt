# Dylan Kenneth Eliot & Google Bard AI

"""

Does similar to the c code in python script but for c++ code. Further extending the languages.

"""


import subprocess
import ctypes

cpp_code = """
#include <iostream>

using namespace std;

void start() {
  cout << "Hello, World!" << endl;
}
"""
proc = subprocess.run(["g++", "-o", "-", "-c"], input=cpp_code.encode("utf-8"), stdout=subprocess.PIPE)
hello = ctypes.CDLL(proc.stdout)
hello.start.restype = ctypes.c_void_p
hello.start()

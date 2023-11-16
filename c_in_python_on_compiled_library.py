# Dylan Kenneth Eliot & Google-Bard-AI

"""
In this file, this basically compiles the code and directly links it instead of worrying about language difference.

But rather compiling to output memory slot of a subprocess, then using that as a means to load the "binary file".
Then do stuff with its insides.

"""

import subprocess
import ctypes

c_code = """
#include <stdio.h>

void start() {
  printf("Hello, World!");
}
"""
proc = subprocess.run(["gcc", "-o", "-", "-c"], input=c_code.encode("utf-8"), stdout=subprocess.PIPE)
hello = ctypes.CDLL(proc.stdout)
hello.start.restype = ctypes.c_void_p
hello.start()

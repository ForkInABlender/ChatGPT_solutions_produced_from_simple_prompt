# Dylan Kenneth Eliot

"""
It took time to find this one and I give that dev its credit as it owns that source code.

It also exists here as it works for [python3.6] to [python3.10].
"""


from cpypp import py_preprocessor
PYPP = py_preprocessor()
PYPP.parse(__file__, __name__)
import sys


"""This source will work in both Python versions"""

#exclude
if len(sys.argv) > 1 and sys.argv[1] == '-d': PYPP.define("debug")
#endexclude

#include "test_0.py"
print(dir())


def main():
   #define f 2
   #expand print(f)
   #ifdef __PYTHON3__
   #expand print("This will work in", "Python 2", "like a charm")
   #else
   print("This will work in Python 3 like a charm")
   #endif

if __name__ == "__main__": main()

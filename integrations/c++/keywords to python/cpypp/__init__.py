"""
PYPP is a preprocessor for language Python versions 2.7+

It support c-style directives in Python with run time and
batch preprocessing, with builtin compiling. You can keep 
your code working in Python 2 and Python 3 forever.

    from cpypp import py_preprocessor
    PYPP = py_preprocessor()
    PYPP.parse(__file__, __name__)

    \"\"\"This source will work in both Python versions\"\"\"

    def main():

       #ifdef __PYTHON2__
       #expand print "This will work in", "Python 2", "like a charm"
       #else
       print("This will work in Python 3 like a charm")
       #endif

    if __name__ == "__main__": main()

$ python2 code_above.py
Only work in Python 2
$ python3 code_above.py
Works from Python 2 and above

Benefits

- Pre-process and run code all in one shot (run mode)
- Batch generate final code files clean of all developer 
  stuff (batch mode) 
- Recompile modules at run time to guarantee version 
  compatibility between Python 2 and Python 3
- Compile final bytecode without temp files
- Produce code that can be run in both Python2 and Python3 
  skipping compatibility problems
- Writing code for debug and has a clean file for production

See https://github.com/wellrats/cpypp for more options
"""
from ._cpypp import py_preprocessor, __version__, __author__, __email__

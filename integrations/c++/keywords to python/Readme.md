# What is this for and why not just compile down?

While compiling down has many benefits, doing it the least number of times as possible is better. In this case, if you don't have to compile at all, all the better.

My reverse ask is "Why compile what will run with a few interchanges to parts layout for language compatibility? Why not use the raw format with the right type castings without compiling?".

# So your goal is to use bits of the higher level abstraction to imitate the lower-level?

Yes. Unicorn-engine plus the following modules will be used for imitation-crabbin' c++ via abstraction, and by implementation. Templating will take place so as to make it possible to fully translatate and move
 c/c++ projects from c/c++ native code to python code with syntactic glue/sugar.

## Why use `unicorn-engine`, `namespaces.py`, `goto.py`, `cpypp` and `ctypes` to imitate c++ in an uncompiled state of execution?

The purpose is to allow for software that would normally be compiled to have a compatibility & translated layer that executes. The idea is that one templates by keyword usage, then refine from their.

The end result should look like c++ code embedded in raw cpython code. This embedded code layer would act as a CPythonC++ version of c++ code to run entirely in cpython or other python like shells. As long as ctypes & _ctypes
 have an implementation.

In additionally would allow for on-the-fly testing of code, leaving errors during compiling behind. As instead you'd pre-execution-wise check everything as it ran in over of each operation, and caring about errors when it
 occurs. And with a set of config-flags, including other needed directories, macros, etc, it is entirely possible to incorporate files used during normal compiling.  

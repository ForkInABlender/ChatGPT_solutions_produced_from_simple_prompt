# What is this for and why not just compile down?

While compiling down has many benefits, doing it the least number of times as possible is better. In this case, if you don't have to compile at all, all the better.

My reverse ask is "Why compile what will run with a few interchanges to parts layout for language compatibility? Why not use the raw format with the right type castings without compiling?".

# So your goal is to use bits of the higher level abstraction to imitate the lower-level?

Yes. Unicorn-engine plus the following modules will be used for imitation-crabbin' c++ via abstraction, and by implementation. Templating will take place so as to make it possible to fully translatate and move
 c/c++ projects from c/c++ native code to python code with syntactic glue/sugar.
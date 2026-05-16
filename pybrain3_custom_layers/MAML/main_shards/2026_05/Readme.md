# What is diskheap.(so/c)?

diskheap.(so/c) is for mapping to and from disk instead of memory.

When an application normally would hit the edge of memory, this is useful as a middle ground process for keeping it all
 running without the application crashing.

# How do i use this at runtime of my binary applications to prevent blowing out the memory table?

``gcc -O2 -shared -fPIC -o diskheap.so diskheap.c -ldl``
<br>then once compiled:<br>
``DISKHEAP_FILE=./diskheap.bin PYTHONMALLOC=malloc LD_PRELOAD=./diskheap.so {{binary go here plus flags}}``<br>
# What does this mean?
That applications that are memory heavy are now disk heavy as an offset to memory.
 Instead of crashing, the app keeps running until end of execution.

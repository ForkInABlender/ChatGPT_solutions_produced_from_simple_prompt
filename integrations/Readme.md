# Are these integrations for python specifically?

Yes. They're the raw form of integration. As long as the types match when used as demonstrated in the integration example .py files.
These are also proof of how to interact with things like kubernetes and its runtime in a dynamic manner. Meaning you can compile bits of how you interact with
 kubernetes.


For example, if you compile the bits needed to get stack information about where your pods are at without having to do a lot of command structuring, you'd use
 bits of its source code to check the same points asked for via command line flags. In addition to this also making it possible to bend a stack with python, 

 pointer dereferencing deeper into the stack can be a problem. For instance, with the script below, one could override the underpinning bits used for others to
  access that "same function". Even though it isn't overriden necessarily on their machine, the backend can be messed with in this mannerism as well.

https://github.com/ForkInABlender/-and-so-sayeth-the-math-nerds-1-1-5-/blob/main/1_add_1_equal_5.py

Not only can numerics be reassigned, but functions entirely swapped from compiled ones to non-compiled versions, etcetera. Meaning if you know what you're doing, you can
easily override a sub function of how kubernetes or docker's backend works.

In addition, one can also use this to dynamically compile portions of dockerd, kubernetes, and docker's source golang for fine grain control and interaction
 with containers they represent and manage.

# What is this for and why not just compile down?

While compiling down has many benefits, doing it the least number of times as possible is better. In this case, if you don't have to compile at all, all the better.

My reverse ask is "Why compile what will run with a few interchanges to parts layout for language compatibility?".

# So your goal is to use bits of the higher level abstraction to imitate the lower-level?

Yes. Unicorn-engine plus the following modules will be used for imitation-crabbin' c++ via abstraction, and by implementation. Templating will take place so as to make it possible to fully translatate and move
 c/c++ projects from c/c++ native code to python code with syntactic glue/sugar.

## Why use `unicorn-engine`, `namespaces.py`, `switch_case.py`, `goto.py`, `pymacros_exec.py`, `cpp_comments_2_func_docs.py` and `ctypes` to imitate c++ in an uncompiled state of execution?

The purpose is to allow for software that would normally be compiled to have a compatibility & translated layer that executes. The idea is that one templates by keyword usage, then refine from their.

The end result should look like c++ code embedded in raw cpython code. This embedded code layer would act as a CPythonC++ version of c++ code to run entirely in cpython or other python like shells. As long as ctypes & _ctypes
 have an implementation.

In additionally would allow for on-the-fly testing of code, leaving errors during compiling behind. As instead you'd pre-execution-wise check everything as it ran in over of each operation, and caring about errors when it
 occurs. And with a set of config-flags, including other needed directories, macros, etc, it is entirely possible to incorporate files used during normal compiling.
This may also mean the results perform 1:1 in terms of execution, utility, and functionality.

Because `namespaces.py`, `switch_case.py`, `goto.py`, & `pymacros_exec.py` enable me to have _ctypes & unicorn-engine operate in an uncompiled state, this allows for even brython.js to use these libraries.

# How would we use `goto`, `label`, `switch`, `case`, `use`, `namespace` in python?

```

from switch import Switch
from goto import label, goto
from namespace import Namespace

x_namespace = Namespace()

@x_namespace
def performAction(choice):
    with Switch(choice) as case:
        while case(1):
            print("Option 1 selected: Hello")
            goto("end_section")
            break
        while case(2):
            print("Option 2 selected: World")
            break
        while case.default:
            print("Invalid choice")
            return
    goto("end_section")


@label
@x_namespace # comment out that namespace annotation, and the function runs as per normal, otherwise it defaults to
def end_section(): # returning a NoneType; which is the equivalent of using `pass` or `nop`.
    print("End of section reached.")

# the above function is now a label that can only be used by namespace 'x_namespace' unless commented out. then it is
#  a function and a 'label' that can be skipped to upon request.


with x_namespace: # using namespace x_namespace
   performAction(int(input("Enter 1 or 2: ")))

```

This works similar to the c++ code normally compiled and then used:

```
#include <iostream>

namespace exampleNamespace {
    void performAction(int choice) {
        switch (choice) {
        case 1:
            std::cout << "Option 1 selected: Hello" << std::endl;
            goto end_section;
        case 2:
            std::cout << "Option 2 selected: World" << std::endl;
            break;
        default:
            std::cout << "Invalid choice" << std::endl;
            return; // Return to avoid hitting the goto label
        }

    end_section:
        std::cout << "End of section reached." << std::endl;
    }
}

int main() {
    int userChoice;
    std::cout << "Enter 1 or 2: ";
    std::cin >> userChoice;
    exampleNamespace::performAction(userChoice);
    return 0;
}
```

This works even for python-{2.7, 3.6, 3.8, 3.9, 3.10} thus far. Meaning it will remain to be forward-compatible, and also capable of supporting backwards compatibility with prior versions of python, etc.

More information on types, go to [my blog](https://brython-fly.blogspot.com/2024/10/when-translating-from-c-to-python-make.html).


# How do c++ comments work in python?

The same way they do in c++, except c++ does not have a `help` function to see the compiled out comments; there are some comments that get filtered out, like the ones around returns, item assignments of value, etcetera. While not perfect, it gets 99% and the remaining .9% to 1% of the work is miniscule or redundant. Do note, though, that all defined comments above, inside, around a function, parameters included, are of that function's `help` comment-tree. 

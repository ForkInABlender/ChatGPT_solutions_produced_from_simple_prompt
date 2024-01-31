# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
Now one can use java, jython, and python code interchangably.

This also allows for making use of java code natively within python code and interoperability. 



https://github.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/blob/2023_09/nasm_to_python_stuff.py
https://github.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/blob/2023_11-1/c%2B%2B_in_python_compiled_library.py
https://github.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/blob/2023_11-1/golang_compiled_library_in_python.py

These three links do similar interaction with dynamic linking steps. Meaning one can also integrate on the fly bits of kubernetes to python.

What even this allows for is operator overriding even of numerical values used as constants by java with the following script 

https://github.com/ForkInABlender/-and-so-sayeth-the-math-nerds-1-1-5-/blob/main/1_add_1_equal_5.py

This also means any value, given the correct type assignment for casting, override even java functions/methods and values used by them.

Please use these wisely.




:)

"""




import jpype
import jpype.imports
from jpype.types import *
jpype.startJVM()

jpype.addClassPath("/home/$USER/Downloads/jython-standalone-2.7.3.jar")

from org.python.util import PythonInterpreter

PythonInterpreter().exec_('print("hi")')

print("Jython loaded successfully!")
jpype.shutdownJVM()

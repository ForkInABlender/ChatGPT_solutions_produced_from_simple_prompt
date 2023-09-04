# Dylan Kenneth Eliot & gpt-3 via openai python module

"""
Objectively, it took one prompt to get the import and include portion. Invariably, it can run commands and follow directions. But, unlike GPT-4 or autogpt, it needs more context awareness and prompt engineering.




"""

import cppyy

cppyy.include("./cpp_code.cpp") # The include I changed from its template/example answer.

# The part I wrote
print(cppyy.gbl.main())

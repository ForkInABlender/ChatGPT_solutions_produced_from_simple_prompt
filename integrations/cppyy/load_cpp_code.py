# Dylan Kenneth Eliot & gpt-3 via openai python module

"""

I reuploaded this one for a different project it is needed for. You may find the cpp_code.cpp in the branch "2023_09" branch on the repository.



"""

import cppyy

cppyy.include("./cpp_code.cpp") # The include I changed from its template/example answer.

# The part I wrote
print(cppyy.gbl.main())

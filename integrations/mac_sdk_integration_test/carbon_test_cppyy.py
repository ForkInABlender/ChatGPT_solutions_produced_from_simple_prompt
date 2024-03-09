# Dylan Kenneth Eliot & GPT-4-Plugins ( Alpha Edition )

"""
This file, with other mac SDK files, allows for "MacTypes" to be imported using their c++ headers.



"""

import cppyy
cppyy.add_include_path("./")
cppyy.cppdef("""
#define TARGET_CPU_X86_64 1
#include "./MacTypes.h"
UInt32 add(UInt32 a, UInt32 b) {
    return a + b;
}
""")
print(f"Result: {cppyy.gbl.add(10, 20)}")

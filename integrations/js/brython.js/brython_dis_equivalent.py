# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This works for disassembling python code within brython.js
 Function reconstruction via this is your own problem to solve.
  Meaning if you retool it some, one could turn that disassembly
   into a useful way to reload bits of cpython. Again, you're on
  your own to reverse engineer such.
"""

import browser
from browser import console
import opcode


def example(a, b):
    return a + b

def simple_disassemble(func):
    code_obj = func.__code__
    bytecode = list(code_obj.co_code)
    instructions = []
    i = 0

    while i < len(bytecode):
        op = bytecode[i]
        op_name = opcode.opname[ord(op)]
        if ord(op) >= opcode.HAVE_ARGUMENT:
            if i + 2 < len(bytecode):
                arg = ord(bytecode[i + 1]) + (ord(bytecode[i + 2]) << 8)
                instructions.append(f"{i:4}: {op_name:<20} {arg}")
                i += 3
            else:
                arg = ord(bytecode[i-3 + 1]) + (ord(bytecode[i-3 + 2]) << 8)
                # If there's not enough bytes left for an argument, break
                instructions.append(f"{i:4}: {op_name:<20} {arg}")
                break
        else:
            instructions.append(f"{i:4}: {op_name}")
            i += 1

    return "\n".join(instructions)

disassembled_code = simple_disassemble(example)
print(disassembled_code)


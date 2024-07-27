# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Adds 2 floats
"""

from unicorn import *
from unicorn.x86_const import *
import struct

mu = Uc(UC_ARCH_X86, UC_MODE_32)
CODE_ADDRESS = 0x1000
DATA_ADDRESS = 0x2000

# Assembly code to add two floats
# fld qword ptr [esp]    ; Load first float from [esp]
# fld qword ptr [esp+8]  ; Load second float from [esp+8]
# faddp st(1), st(0)     ; Add st(1) and st(0), store result in st(1)
# fstp qword ptr [esp]   ; Store result back to [esp]
X86_CODE32 = b'\xdd\x04\x24\xdd\x44\x24\x08\xde\xc1\xdd\x1c\x24'

mu.mem_map(CODE_ADDRESS, 4096)
mu.mem_map(DATA_ADDRESS, 4096)
mu.mem_write(CODE_ADDRESS, X86_CODE32)

def emulate_and_get_result(a, b):
    FLOAT_1 = struct.pack('d', a)
    FLOAT_2 = struct.pack('d', b)
    mu.mem_write(DATA_ADDRESS, FLOAT_1 + FLOAT_2)
    mu.reg_write(UC_X86_REG_ESP, DATA_ADDRESS)
    try:
        mu.emu_start(CODE_ADDRESS, CODE_ADDRESS + len(X86_CODE32))
        result = mu.mem_read(DATA_ADDRESS, 8)
        result = struct.unpack('d', result)[0]
    except UcError as e:
        print(f"Unicorn Engine Error: {e}")
    return result

print(emulate_and_get_result(3.5, 2.5))

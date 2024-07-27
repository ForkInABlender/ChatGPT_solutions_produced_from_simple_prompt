# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Adds 2 integers
"""

from unicorn import *
from unicorn.x86_const import *

# Memory address where emulation starts
ADDRESS = 0x0001000
STACK_ADDR = 0x0002000
STACK_SIZE = 4096

X86_CODE32 = b"\x8b\x44\x24\x04"  # mov eax, [esp + 4]
X86_CODE32 += b"\x03\x44\x24\x08"  # add eax, [esp + 8]
X86_CODE32 += b"\xc3"              # ret

mu = Uc(UC_ARCH_X86, UC_MODE_32)
mu.mem_map(ADDRESS, STACK_SIZE)
mu.mem_map(STACK_ADDR, STACK_SIZE)
mu.mem_write(ADDRESS, X86_CODE32)

def emulate_add(a, b):
    stack_base = STACK_ADDR + STACK_SIZE - 16
    #print(f"Writing value {a} to address {hex(stack_base + 4)}")
    #print(f"Writing value {b} to address {hex(stack_base + 8)}")
    mu.mem_write(stack_base + 4, a.to_bytes(4, byteorder='little'))  # a
    mu.mem_write(stack_base + 8, b.to_bytes(4, byteorder='little'))  # b
    mu.reg_write(UC_X86_REG_ESP, stack_base)
    #print(f"ESP set to {hex(mu.reg_read(UC_X86_REG_ESP))}:")
    mu.reg_write(UC_X86_REG_EIP, ADDRESS)
    try:
        mu.emu_start(ADDRESS, ADDRESS + len(X86_CODE32))
    except UcError as e:
        pass
    try:
        esp_plus_4 = mu.mem_read(stack_base + 4, 4)
        esp_plus_8 = mu.mem_read(stack_base + 8, 4)
        #print(f"Value at ESP+4: {int.from_bytes(esp_plus_4, 'little')}")
        #print(f"Value at ESP+8: {int.from_bytes(esp_plus_8, 'little')}")
    except UcError as e:
        print(f"Memory read error: {e}")
    result = mu.reg_read(UC_X86_REG_EAX)
    return result

print("Result: {}".format(emulate_add(2, 3)))

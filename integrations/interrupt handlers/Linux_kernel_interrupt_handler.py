#Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""

This is the linix interrupt handler.

This version works even on android through UserLAnd apk running ubuntu jammy.
The windows kernel interrupt handler can be equally written this way if you're aiming to
run dockerd and docker on Android but without loading an entire kernel over your existing one.

The test file t2.py and terminal output have been also uploaded for reference.
"""


from unicorn import *
from unicorn.x86_const import *

# x86 machine code to simulate an interrupt
# xor eax, eax
# mov eax, 1
# add eax, 3
# int 0x80  ; Interrupt instruction for syscall
X86_CODE32 = b"\x31\xc0" \
             b"\xb8\x01\x00\x00\x00" \
             b"\x05\x03\x00\x00\x00" \
             b"\xcd\x80"

# memory address where emulation starts
ADDRESS = 0x1000000

# Initialize emulator in X86-32bit mode
mu = Uc(UC_ARCH_X86, UC_MODE_32)

# map 2MB memory for this emulation
mu.mem_map(ADDRESS, 2 * 1024 * 1024)

# write machine code to be emulated to memory
mu.mem_write(ADDRESS, X86_CODE32)

# initialize machine registers
mu.reg_write(UC_X86_REG_EAX, 0x0)

# Hook for interrupt
def hook_intr(uc, intno, user_data):
    if intno == 0x80:  # Software interrupt for syscalls
        print("Syscall Interrupt Triggered")
        eax = uc.reg_read(UC_X86_REG_EAX)
        print(f"EAX before interrupt: {eax}")
        # You can modify register or memory based on your needs here
        uc.reg_write(UC_X86_REG_EAX, eax + 5)  # example modification

mu.hook_add(UC_HOOK_INTR, hook_intr)

# emulate code in infinite time & unlimited instructions
mu.emu_start(ADDRESS, ADDRESS + len(X86_CODE32))

# now print out the EAX register
print("EAX after interrupt = 0x%x" % mu.reg_read(UC_X86_REG_EAX))

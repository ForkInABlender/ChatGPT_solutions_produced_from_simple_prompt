import ctypes
import unicorn

# define architecture for x86_64
ARCH_X86_64 = unicorn.UC_ARCH_X86
MODE_64 = unicorn.UC_MODE_64

# initialize emulator in X86_64 mode
emulator = unicorn.Uc(ARCH_X86_64, MODE_64)

# write machine code to be emulated
machine_code = b"\x41\x4a" # inc edx, dec edx in x86_64

# map 2 MB of memory for this emulation
emulator.mem_map(0x1000, 2 * 1024 * 1024)

# write machine code to memory
emulator.mem_write(0x1000, machine_code)

# set the PC to start execution from
emulator.reg_write(unicorn.x86_const.UC_X86_REG_RIP, 0x1000)

# Callback for handling interrupt request
def hook_intr(uc, intno, data):
    if intno == 0x2E:
        # Handle Windows system call
        print("Handling Windows system call")
        return True
    return False

# Register the interrupt callback
emulator.hook_add(unicorn.UC_HOOK_INTR, hook_intr)

# start emulation, code will be executed for 1000 instructions
emulator.emu_start(0x1000, 0x100 + len(machine_code), 0, 1000)

# read register value after emulation
result = emulator.reg_read(unicorn.x86_const.UC_X86_REG_EDX)

print("Result after execution: %x" % result)

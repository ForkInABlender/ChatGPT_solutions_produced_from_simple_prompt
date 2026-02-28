import sys
import struct
from unicorn import *
from unicorn.x86_const import *

# -----------------------------
# Constants
# -----------------------------

SETUP_ADDR = 0x90000
KERNEL_ADDR = 0x100000
MEM_SIZE = 128 * 1024 * 1024
BOOT_HEADER_OFFSET = 0x1f1

# -----------------------------
# Parse bzImage
# -----------------------------

def parse_bzimage(path):
    with open(path, "rb") as f:
        data = f.read()

    if data[0x202:0x206] != b'HdrS':
        raise ValueError("Invalid bzImage (missing HdrS)")

    setup_sects = data[BOOT_HEADER_OFFSET]
    if setup_sects == 0:
        setup_sects = 4

    setup_size = (setup_sects + 1) * 512
    setup = data[:setup_size]
    kernel = data[setup_size:]

    return setup, kernel

# -----------------------------
# Interrupt Hook
# -----------------------------

def hook_intr(uc, intno, user_data):
    print(f"[INT] 0x{intno:02x}")

    # Minimal fake BIOS responses
    if intno == 0x15:
        # fake success
        uc.reg_write(UC_X86_REG_CF, 0)
        return

    if intno == 0x10:
        # video BIOS
        uc.reg_write(UC_X86_REG_CF, 0)
        return

    raise Exception(f"Unhandled interrupt 0x{intno:x}")

# -----------------------------
# CR0 Hook (detect protected mode)
# -----------------------------

def hook_code(uc, address, size, user_data):
    # Detect CR0 writes
    try:
        opcode = uc.mem_read(address, size)
    except:
        return

    # naive detection for "mov cr0"
    if opcode.startswith(b"\x0f\x22"):  # MOV to CRx
        print("[*] Possible CR write detected")

# -----------------------------
# Transition to 32-bit mode
# -----------------------------

def switch_to_protected(uc16):
    print("[*] Switching emulator to 32-bit mode")

    # Save state
    regs = {}
    for reg in [
        UC_X86_REG_EAX, UC_X86_REG_EBX,
        UC_X86_REG_ECX, UC_X86_REG_EDX,
        UC_X86_REG_ESI, UC_X86_REG_EDI,
        UC_X86_REG_EBP, UC_X86_REG_ESP,
        UC_X86_REG_EIP
    ]:
        regs[reg] = uc16.reg_read(reg)

    mem = uc16.mem_read(0, MEM_SIZE)

    # Create new emulator in 32-bit mode
    uc32 = Uc(UC_ARCH_X86, UC_MODE_32)
    uc32.mem_map(0, MEM_SIZE)
    uc32.mem_write(0, mem)

    for reg, val in regs.items():
        uc32.reg_write(reg, val)

    return uc32

# -----------------------------
# Main
# -----------------------------

def main(path):
    setup, kernel = parse_bzimage(path)

    uc = Uc(UC_ARCH_X86, UC_MODE_16)
    uc.mem_map(0, MEM_SIZE)

    uc.mem_write(SETUP_ADDR, setup)
    uc.mem_write(KERNEL_ADDR, kernel)

    # Real mode register setup
    uc.reg_write(UC_X86_REG_CS, 0x9000)
    uc.reg_write(UC_X86_REG_IP, 0x0000)

    uc.reg_write(UC_X86_REG_DS, 0x9000)
    uc.reg_write(UC_X86_REG_ES, 0x9000)
    uc.reg_write(UC_X86_REG_SS, 0x9000)
    uc.reg_write(UC_X86_REG_SP, 0x8000)

    uc.hook_add(UC_HOOK_INTR, hook_intr)
    uc.hook_add(UC_HOOK_CODE, hook_code)

    print("[*] Starting 16-bit execution")

    try:
        uc.emu_start(SETUP_ADDR, SETUP_ADDR + len(setup))
    except UcError as e:
        print("[!] Emulation stopped:", e)

    print("[*] Done")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python loader.py /path/to/bzImage")
        sys.exit(1)

    main(sys.argv[1])

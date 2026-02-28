import ctypes
import struct
import sys
import os

# Constants
BZIMAGE_LOAD_ADDR = 0x100000   # 1 MB
SETUP_ADDR = 0x90000           # Real-mode setup code loaded here
BOOT_HEADER_OFFSET = 0x1f1     # Offset of setup_sects in the boot header
MEMORY_SIZE = 64 * 1024 * 1024 # 64 MB simulated RAM

# ------------------ BootParams Struct ------------------

class BootParams(ctypes.LittleEndianStructure):
    _fields_ = [
        ("setup_sects", ctypes.c_ubyte),
        ("root_flags", ctypes.c_uint16),
        ("sys_size", ctypes.c_uint16),
        ("ram_size", ctypes.c_uint16),
        ("vid_mode", ctypes.c_uint16),
        ("root_dev", ctypes.c_uint16),
        ("boot_flag", ctypes.c_uint16),
        ("jump", ctypes.c_uint16),
        ("header", ctypes.c_uint32),             # == 0x53726448 ("HdrS")
        ("version", ctypes.c_uint16),
        ("realmode_swtch", ctypes.c_uint32),
        ("start_sys_seg", ctypes.c_uint16),
        ("kernel_version", ctypes.c_uint16),
        ("type_of_loader", ctypes.c_uint8),
        ("loadflags", ctypes.c_uint8),
        ("setup_move_size", ctypes.c_uint16),
        ("code32_start", ctypes.c_uint32),
        ("ramdisk_image", ctypes.c_uint32),
        ("ramdisk_size", ctypes.c_uint32),
        ("bootsect_kludge", ctypes.c_uint32),
        ("heap_end_ptr", ctypes.c_uint16),
        ("ext_loader_ver", ctypes.c_uint8),
        ("ext_loader_type", ctypes.c_uint8),
        ("cmd_line_ptr", ctypes.c_uint32),
        ("initrd_addr_max", ctypes.c_uint32),
        ("kernel_alignment", ctypes.c_uint32),
        ("relocatable_kernel", ctypes.c_uint8),
        ("min_alignment", ctypes.c_uint8),
        ("xloadflags", ctypes.c_uint16),
        ("cmdline_size", ctypes.c_uint32),
        ("hardware_subarch", ctypes.c_uint32),
        ("hardware_subarch_data", ctypes.c_uint64),
        ("payload_offset", ctypes.c_uint32),
        ("payload_length", ctypes.c_uint32),
        ("setup_data", ctypes.c_uint64),
        ("pref_address", ctypes.c_uint64),
        ("init_size", ctypes.c_uint32),
        ("handover_offset", ctypes.c_uint32),
        ("kernel_info_offset", ctypes.c_uint32),
    ]

# ------------------ Simulated Memory ------------------

class SimulatedMemory:
    def __init__(self, size: int = MEMORY_SIZE):
        self.size = size
        self.mem = (ctypes.c_ubyte * size)()
        self.base_addr = ctypes.addressof(self.mem)

    def load_bytes(self, addr: int, data: bytes):
        if addr + len(data) > self.size:
            raise MemoryError("Attempt to load beyond memory limits")
        ctypes.memmove(self.base_addr + addr, data, len(data))

    def read_bytes(self, addr: int, size: int) -> bytes:
        if addr + size > self.size:
            raise MemoryError("Attempt to read beyond memory limits")
        return bytes((self.mem)[addr:addr + size])

    def read_struct(self, addr: int, struct_type):
        data = self.read_bytes(addr, ctypes.sizeof(struct_type))
        return struct_type.from_buffer_copy(data)

# ------------------ bzImage Loader Logic ------------------

def parse_bzimage(image_bytes: bytes):
    if len(image_bytes) < 0x206:
        raise ValueError("bzImage too small to contain a valid header")

    hdr_signature = image_bytes[0x202:0x206]
    if hdr_signature != b'HdrS':
        raise ValueError("bzImage does not contain a valid boot header (missing HdrS signature)")

    setup_sects = image_bytes[BOOT_HEADER_OFFSET]
    if setup_sects == 0:
        setup_sects = 4  # Default per boot protocol

    setup_bytes = (setup_sects + 1) * 512
    if len(image_bytes) < setup_bytes:
        raise ValueError("bzImage is truncated: setup size exceeds file size")

    setup_data = image_bytes[:setup_bytes]
    kernel_data = image_bytes[setup_bytes:]

    return setup_data, kernel_data

def prepare_bzimage_memory(image_path: str) -> SimulatedMemory:
    if not os.path.isfile(image_path):
        raise FileNotFoundError(f"bzImage file not found: {image_path}")

    with open(image_path, 'rb') as f:
        bzimage_data = f.read()

    setup_data, kernel_data = parse_bzimage(bzimage_data)

    memory = SimulatedMemory()

    # Load setup and kernel into simulated memory
    memory.load_bytes(SETUP_ADDR, setup_data)
    memory.load_bytes(BZIMAGE_LOAD_ADDR, kernel_data)

    print(f"[+] Setup code loaded at 0x{SETUP_ADDR:05x} ({len(setup_data)} bytes)")
    print(f"[+] Kernel loaded at 0x{BZIMAGE_LOAD_ADDR:06x} ({len(kernel_data)} bytes)")

    # Read the BootParams from the start of the setup area (not from offset!)
    boot_params = memory.read_struct(SETUP_ADDR, BootParams)

    if boot_params.header != 0x53726448:
        raise ValueError("Invalid boot_params header signature")

    print("\n[=] BootParams extracted:")
    print(f"    setup_sects        : {boot_params.setup_sects}")
    print(f"    version            : 0x{boot_params.version:04x}")
    print(f"    code32_start       : 0x{boot_params.code32_start:08x}")
    print(f"    ramdisk_image      : 0x{boot_params.ramdisk_image:08x}")
    print(f"    ramdisk_size       : {boot_params.ramdisk_size}")
    print(f"    type_of_loader     : {boot_params.type_of_loader}")
    print(f"    cmd_line_ptr       : 0x{boot_params.cmd_line_ptr:08x}")
    print(f"    cmdline_size       : {boot_params.cmdline_size}")
    print(f"    kernel_alignment   : {boot_params.kernel_alignment}")
    print(f"    relocatable_kernel : {boot_params.relocatable_kernel}")
    print(f"    xloadflags         : 0x{boot_params.xloadflags:04x}")
    print(f"    pref_address       : 0x{boot_params.pref_address:016x}")
    print(f"    init_size          : 0x{boot_params.init_size:08x}")
    print(f"    handover_offset    : 0x{boot_params.handover_offset:08x}")
    print(f"    kernel_info_offset : 0x{boot_params.kernel_info_offset:08x}")

    return memory

# ------------------ Main Entrypoint ------------------

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} /path/to/bzImage")
        sys.exit(1)

    image_path = sys.argv[1]

    try:
        memory = prepare_bzimage_memory(image_path)
        print("\n[✓] bzImage memory preparation complete.")
    except Exception as ex:
        print(f"[!] Error: {ex}")
        sys.exit(1)

if __name__ == '__main__':
    main()

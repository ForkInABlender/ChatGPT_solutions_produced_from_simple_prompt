# Dylan Kenneth Eliot

"""
This will later be used with unicorn engine to run windows applications on linux without the need for wine.

later on, similar will be done for linux binaries so they work with windows without WSL, hypervisor, or other virtualization "technology". 

"""

import ctypes
from enum import IntEnum

class COFFHeader(ctypes.Structure):
		_fields_ = [
				("magic", ctypes.c_char * 4),
				("machine", ctypes.c_uint16),
				("numberOfSections", ctypes.c_uint16),
				("timeDateStamp", ctypes.c_uint32),
				("pointerToSymbolTable", ctypes.c_uint32),
				("numberOfSymbols", ctypes.c_uint32),
				("sizeOfOptionalHeader", ctypes.c_uint16),
				("characteristics", ctypes.c_uint16)
		]

# Define the pe_subsystem enumeration
class PeSubsystem(IntEnum):
    IMAGE_SUBSYSTEM_UNKNOWN = 0x0
    IMAGE_SUBSYSTEM_NATIVE = 0x1
    IMAGE_SUBSYSTEM_WINDOWS_GUI = 0x2
    IMAGE_SUBSYSTEM_WINDOWS_CUI = 0x3
    IMAGE_SUBSYSTEM_OS2_CUI = 0x5
    IMAGE_SUBSYSTEM_POSIX_CUI = 0x7
    IMAGE_SUBSYSTEM_NATIVE_WINDOWS = 0x8
    IMAGE_SUBSYSTEM_WINDOWS_CE_GUI = 0x9
    IMAGE_SUBSYSTEM_EFI_APPLICATION = 0xa
    IMAGE_SUBSYSTEM_EFI_BOOT_SERVICE_DRIVER = 0xb
    IMAGE_SUBSYSTEM_EFI_RUNTIME_DRIVER = 0xc
    IMAGE_SUBSYSTEM_EFI_ROM = 0xd
    IMAGE_SUBSYSTEM_XBOX = 0xe
    IMAGE_SUBSYSTEM_WINDOWS_BOOT_APPLICATION = 0x10

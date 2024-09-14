# Dylan Kenneth Eliot

"""


"""

import ctypes

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

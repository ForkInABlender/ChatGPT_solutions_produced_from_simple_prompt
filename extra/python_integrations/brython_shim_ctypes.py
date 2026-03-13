from __future__ import annotations
# Dylan Kenneth Eliot - Brython Optimized _ctypes shim

import array as _array
import struct as _struct
import sys as _sys
import os as _os
import typing as _typing
import weakref
from collections import OrderedDict
from browser import window

# ---------------------------------------------------------------------------
# Platform constants for Brython (typically 32-bit JS execution context)
# ---------------------------------------------------------------------------
_PTR_SZ   = 4  # Standard for Brython/JS
_PTR_FMT  = 'I'
_ENDIAN   = '<' 

__version__             = '1.1.0'
CTYPES_MAX_ARGCOUNT     = 1024
FUNCFLAG_CDECL          = 0x1
FUNCFLAG_PYTHONAPI      = 0x4
FUNCFLAG_USE_ERRNO      = 0x8
FUNCFLAG_USE_LASTERROR  = 0x10
RTLD_LOCAL              = 0
RTLD_GLOBAL             = 0x100
SIZEOF_TIME_T           = 4

_pointer_type_cache: dict = {}

class ArgumentError(Exception): pass

# ---------------------------------------------------------------------------
# Brython-Specific Type Table Build
# ---------------------------------------------------------------------------
_TYPE_TABLE: dict[str, tuple[str, int, int]] = {}

def _build_type_table():
    _StructError = getattr(_struct, 'error', Exception)
    # Define scalars with fallbacks for Brython's VFS _struct
    scalars = [('b','b'),('B','B'),('h','h'),('H','H'),('i','i'),('I','I'),
               ('q','q'),('Q','Q'),('f','f'),('d','d'),('?','?'),('c','c')]
    for code, fmt in scalars:
        try:
            sz = _struct.calcsize('=' + fmt)
            aln = sz
            _TYPE_TABLE[code] = (fmt, sz, aln)
        except:
            _TYPE_TABLE[code] = ('B', 1, 1) # Minimalist fallback

    # Pointer-width types
    for code in ('P', 'z', 'Z', 'O'):
        _TYPE_TABLE[code] = (_PTR_FMT, _PTR_SZ, _PTR_SZ)
    _TYPE_TABLE['u'] = ('H', 2, 2)
    _TYPE_TABLE['g'] = ('d', 8, 8)

_build_type_table()

def _type_info(code: str):
    return _TYPE_TABLE.get(code, ('B', 1, 1))

# ---------------------------------------------------------------------------
# The Memory Engine: Brython Virtual Addressing
# ---------------------------------------------------------------------------

class _RawBuf:
    __slots__ = ('_arr', '_size', 'address')

    def __init__(self, size: int, init: bytes = None):
        if size == 0: size = 1
        self._arr  = _array.array('B', [0]*size)
        self._size = size
        # BRYTHON FIX: Use id() as the virtual memory address
        self.address = id(self._arr)

        if init:
            n = min(len(init), size)
            for i in range(n):
                self._arr[i] = init[i]

    def read(self, offset: int = 0, size: int | None = None) -> bytes:
        end = self._size if size is None else offset + size
        return bytes(self._arr[offset:end])

    def write(self, data, offset: int = 0):
        if isinstance(data, (bytes, bytearray)):
            for i, b in enumerate(data):
                if offset + i < self._size:
                    self._arr[offset + i] = b

    def pack_into(self, fmt, offset, *values):
        full_fmt = _ENDIAN + fmt
        sz = _struct.calcsize(full_fmt)
        # Ensure buffer safety for Brython's struct implementation
        if offset + sz > len(self._arr):
            self._arr.extend([0] * (offset + sz - len(self._arr)))
            self._size = len(self._arr)
        _struct.pack_into(full_fmt, self._arr, offset, *values)

    def unpack_from(self, fmt: str, offset: int):
        return _struct.unpack_from(_ENDIAN + fmt, self._arr, offset)

    def resize(self, new_size: int):
        if new_size <= self._size: return
        self._arr.extend([0] * (new_size - self._size))
        self._size = new_size

    def __len__(self) -> int: return self._size

_addr_to_buf: dict[int, _RawBuf] = {}

def _alloc(size: int, init: bytes | None = None) -> _RawBuf:
    buf = _RawBuf(size, init)
    _addr_to_buf[buf.address] = buf
    return buf

# ---------------------------------------------------------------------------
# Required Functions (The "Must Have" List from _ctypes_pure.py)
# ---------------------------------------------------------------------------

def sizeof(obj) -> int:
    tp = obj if isinstance(obj, type) else type(obj)
    return getattr(tp, '_size_', 0)

def alignment(obj) -> int:
    tp = obj if isinstance(obj, type) else type(obj)
    return getattr(tp, '_align_', 1)

def addressof(obj) -> int:
    if hasattr(obj, '_buffer_'):
        return obj._buffer_.address + getattr(obj, '_offset_', 0)
    raise TypeError('invalid type for addressof')

def byref(obj, offset: int = 0):
    return _CArgObject(obj, offset)

class _CArgObject:
    def __init__(self, obj, offset: int = 0):
        self._obj = obj
        self._offset = offset

# Function Pointer and Library Logic
_vcallables = {}
_vaddr_counter = 0x7FFF0000

def _next_vaddr():
    global _vaddr_counter
    _vaddr_counter += 8
    return _vaddr_counter

def register_library(name: str, symbols: dict[str, callable]) -> int:
    h = id(name)
    for sym, func in symbols.items():
        addr = _next_vaddr()
        _vcallables[addr] = func
        # Map to JS window for Brython visibility
        setattr(window, f"__vlib_{sym}", func)
    return h

# ---------------------------------------------------------------------------
# Base Classes & Data Types (Re-Implementation for Brython Context)
# ---------------------------------------------------------------------------

class _CData:
    _b_needsfree_ = False
    _buffer_ = None
    _offset_ = 0

class _SimpleCData(_CData, metaclass=PyCSimpleType):
    def __init__(self, value=None):
        self._buffer_ = _alloc(self._size_)
        if value is not None: self.value = value
    
    @property
    def value(self):
        # Implementation of _read_value from your file
        return self._buffer_.unpack_from(self._fmt_, self._offset_)[0]
    
    @value.setter
    def value(self, v):
        self._buffer_.pack_into(self._fmt_, self._offset_, v)

# ... [Include Pointer, Structure, Union, and Array logic as per your file] ...

# ---------------------------------------------------------------------------
# External Linkage Addresses (Expected by your integrations)
# ---------------------------------------------------------------------------
_cast_addr       = _next_vaddr()
_memmove_addr    = _next_vaddr()
_memset_addr     = _next_vaddr()
_string_at_addr  = _next_vaddr()
_wstring_at_addr = _next_vaddr()

def _impl_cast(obj, obj2, typ): return cast(obj2 if obj2 else obj, typ)
_vcallables[_cast_addr] = _impl_cast

# ---------------------------------------------------------------------------
# Final Export Surface
# ---------------------------------------------------------------------------
__all__ = [
    'ArgumentError', 'Array', 'CFuncPtr', 'CTYPES_MAX_ARGCOUNT',
    'FUNCFLAG_CDECL', 'FUNCFLAG_PYTHONAPI', 'FUNCFLAG_USE_ERRNO',
    'FUNCFLAG_USE_LASTERROR', 'POINTER', 'PyObj_FromPtr', 'Py_DECREF',
    'Py_INCREF', 'RTLD_GLOBAL', 'RTLD_LOCAL', 'SIZEOF_TIME_T',
    'Structure', 'Union', '_Pointer', '_SimpleCData', '_cast_addr',
    '_memmove_addr', '_memset_addr', '_pointer_type_cache',
    '_string_at_addr', '_unpickle', '_wstring_at_addr',
    'addressof', 'alignment', 'buffer_info', 'byref',
    'call_cdeclfunction', 'call_function', 'dlclose', 'dlopen', 'dlsym',
    'get_errno', 'pointer', 'resize', 'set_errno', 'sizeof',
    'register_library',
]

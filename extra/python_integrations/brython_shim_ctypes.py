from __future__ import annotations
# Dylan Kenneth Eliot - Brython Optimized _ctypes shim

"""
This is the brython shim being tested on premise that it will work during testing & has worked for cpython similarly.

"""

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

# ---------------------------------------------------------------------------
# PyCSimpleType (Metaclass must come FIRST)
# ---------------------------------------------------------------------------

class PyCSimpleType(type):
    def __new__(mcs, name, bases, ns):
        cls  = super().__new__(mcs, name, bases, ns)
        code = ns.get('_type_')
        if code is not None:
            fmt, sz, aln = _type_info(code)
            cls._size_   = sz
            cls._align_  = aln
            cls._fmt_    = fmt
        return cls

    def __mul__(cls, count: int):
        return _make_array_type(cls, count)

    def from_address(cls, address: int):
        buf = _buf_containing_addr(address)
        if buf is None:
            buf = _ExternalBuf(address, cls._size_)
        inst = cls.__new__(cls)
        inst._buffer_       = buf
        inst._offset_       = address - buf.address
        inst._b_needsfree_  = False
        return inst

    def from_buffer(cls, source):
        # Brython handles byte-like objects slightly differently; 
        # ensure we get a real byte array
        ba   = bytes(source)
        buf  = _alloc(len(ba), ba)
        inst = cls.__new__(cls)
        inst._buffer_      = buf
        inst._offset_      = 0
        inst._b_needsfree_ = True
        return inst

    def from_buffer_copy(cls, source):
        data = bytes(source)[:cls._size_]
        inst = cls()
        inst._buffer_.write(data.ljust(cls._size_, b'\x00'))
        return inst

    def in_dll(cls, lib, name: str):
        return cls.from_address(_vlib_sym(lib, name))

    def from_param(cls, value):
        if isinstance(value, cls):
            return value
        if hasattr(cls, '_type_') and cls._type_ in ('z', 'Z', 'P') and value is None:
            return cls()
        try:
            return cls(value)
        except (TypeError, ValueError) as e:
            raise TypeError(f"wrong type: expected {cls.__name__}") from e

# ---------------------------------------------------------------------------
# _SimpleCData (Now defined with PyCSimpleType available)
# ---------------------------------------------------------------------------

class _SimpleCData(_CData, metaclass=PyCSimpleType):
    def __init__(self, value=None):
        self._buffer_      = _alloc(self._size_)
        self._offset_      = 0
        self._b_needsfree_ = True
        if value is not None:
            self.value = value

    @property
    def value(self):
        code = self._type_
        fmt  = self._fmt_
        off  = self._offset_

        # Brython/JS logic for specific C-types
        if code == 'c':
            return self._buffer_.read(off, 1)
        if code == 'u':
            wsz = self._size_
            raw = self._buffer_.read(off, wsz)
            enc = 'utf-16-le' if wsz == 2 else 'utf-32-le'
            return raw.decode(enc)
        if code == 'z':
            addr = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            return _read_cstring(addr) if addr != 0 else None
        if code == 'O':
            idx = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            return _py_object_store.get(idx)

        return self._buffer_.unpack_from(fmt, off)[0]

    @value.setter
    def value(self, v):
        code = self._type_
        fmt  = self._fmt_
        off  = self._offset_
        
        # Handle string pointers by allocating temporary virtual buffers
        if code in ('z', 'Z', 'P'):
            if v is None:
                self._buffer_.pack_into(_PTR_FMT, off, 0)
                return
            if isinstance(v, (bytes, str)):
                b = v.encode() + b'\x00' if isinstance(v, str) else v + b'\x00'
                tmp = _alloc(len(b), b)
                self._buffer_.pack_into(_PTR_FMT, off, tmp.address)
                return

        try:
            self._buffer_.pack_into(fmt, off, v)
        except Exception as e:
            raise OverflowError(str(e))


class PyCStructType(type):

    def __new__(mcs, name, bases, ns):
        cls = super().__new__(mcs, name, bases, ns)
        if '_fields_' in ns:
            PyCStructType._apply_fields(cls)
        return cls

    @staticmethod
    def _apply_fields(cls):
        descs, sz, aln = _layout_fields(
            cls._fields_,
            pack=getattr(cls, '_pack_', None),
            is_union=False)
        cls._size_  = sz
        cls._align_ = aln
        for name, desc in descs.items():
            setattr(cls, name, desc)

    def __setattr__(cls, name, value):
        super().__setattr__(name, value)
        if name == '_fields_':
            PyCStructType._apply_fields(cls)

    def __mul__(cls, count: int):
        return _make_array_type(cls, count)

    def from_address(cls, address: int):
        buf  = _buf_containing_addr(address) or _ExternalBuf(address, cls._size_)
        inst = cls.__new__(cls)
        inst._buffer_      = buf
        inst._offset_      = address - buf.address
        inst._b_needsfree_ = False
        return inst

    def from_buffer(cls, source):
        ba  = bytes(source)
        buf = _alloc(len(ba), ba)
        inst = cls.__new__(cls)
        inst._buffer_      = buf
        inst._offset_      = 0
        inst._b_needsfree_ = True
        return inst

class Structure(_CData, metaclass=PyCStructType):
    """Base class for ctypes Structure types."""

    _size_  = 0
    _align_ = 1

    def __init__(self, *args, **kwargs):
        self._buffer_      = _alloc(max(self._size_, 1))
        self._offset_      = 0
        self._b_needsfree_ = True
        fields = getattr(self, '_fields_', [])
        for spec, value in zip(fields, args):
            setattr(self, spec[0], value)
        for fname, value in kwargs.items():
            setattr(self, fname, value)

    def __repr__(self):
        return f'<{type(self).__name__} at 0x{addressof(self):016x}>'

class UnionType(PyCStructType):

    def __new__(mcs, name, bases, ns):
        cls = type.__new__(mcs, name, bases, ns)
        if '_fields_' in ns:
            UnionType._apply_fields(cls)
        return cls

    @staticmethod
    def _apply_fields(cls):
        descs, sz, aln = _layout_fields(
            cls._fields_,
            pack=getattr(cls, '_pack_', None),
            is_union=True)
        cls._size_  = sz
        cls._align_ = aln
        for name, desc in descs.items():
            setattr(cls, name, desc)

    def __setattr__(cls, name, value):
        type.__setattr__(cls, name, value)
        if name == '_fields_':
            UnionType._apply_fields(cls)


class Union(_CData, metaclass=UnionType):
    """Base class for ctypes Union types."""

    _size_  = 0
    _align_ = 1

    def __init__(self, *args, **kwargs):
        self._buffer_      = _alloc(max(self._size_, 1))
        self._offset_      = 0
        self._b_needsfree_ = True
        fields = getattr(self, '_fields_', [])
        if args and fields:
            setattr(self, fields[0][0], args[0])
        for fname, value in kwargs.items():
            setattr(self, fname, value)

    def __repr__(self):
        return f'<{type(self).__name__} at 0x{addressof(self):016x}>'

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

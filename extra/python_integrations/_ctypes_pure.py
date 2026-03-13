# Dylan Kenneth Eliot

"""
This is a _ctypes module shim.

This is useful for when compiling _ctypes down just ain't gone cut it.

"""

"""
_ctypes_pure.py
===============
Pure-Python re-implementation of CPython's _ctypes C extension.

Dependencies: ONLY Python standard library
  - array   : array.array('B') gives a real OS-level buffer address via
              .buffer_info() and works directly with struct.pack_into /
              struct.unpack_from.
  - struct  : all field packing / unpacking / size / alignment.
  - mmap    : (optional) large anonymous allocations.
  - sys, os, weakref, collections : bookkeeping.

NO ctypes, NO cffi, NO C extensions of any kind.

FFI note
--------
dlopen / dlsym / call_function require an OS-level bridge that only exists
in compiled C code.  This module provides a *virtual library registry*:
Python callables can be registered under a (libname, symbol) key and called
through the normal CFuncPtr machinery.  Code that only uses the Python type
system (Structure, Union, Array, _Pointer, _SimpleCData, sizeof, addressof,
byref, …) works fully; code that needs to call real native libraries must
supply Python-callable wrappers via register_library().
"""
from __future__ import annotations

import array  as _array
import struct as _struct
import sys    as _sys
import os     as _os
import typing as _typing

import weakref
from collections import OrderedDict

# ---------------------------------------------------------------------------
# Platform constants
# ---------------------------------------------------------------------------
_PTR_SZ   = _struct.calcsize('P')          # 8 on 64-bit, 4 on 32-bit
_PTR_FMT  = 'Q' if _PTR_SZ == 8 else 'I'  # unsigned int wide enough for a pointer
_ENDIAN   = '<' if _sys.byteorder == 'little' else '>'

# Module-level constants (match real _ctypes)
__version__             = '1.1.0'
CTYPES_MAX_ARGCOUNT     = 1024
FUNCFLAG_CDECL          = 0x1
FUNCFLAG_PYTHONAPI      = 0x4
FUNCFLAG_USE_ERRNO      = 0x8
FUNCFLAG_USE_LASTERROR  = 0x10
RTLD_LOCAL              = 0
RTLD_GLOBAL             = 0x100
SIZEOF_TIME_T           = 8 if _PTR_SZ == 8 else 4

_pointer_type_cache: dict = {}

# ---------------------------------------------------------------------------
# ArgumentError
# ---------------------------------------------------------------------------
class ArgumentError(Exception):
    pass

# ---------------------------------------------------------------------------
# Type-code table
# Maps _type_ character → (struct_fmt, size, alignment)
# We derive sizes & alignments purely from the struct module.
# ---------------------------------------------------------------------------
_TYPE_TABLE: dict[str, tuple[str, int, int]] = {}

def _build_type_table() -> None:
    # Standard scalars: struct format determines size & natural alignment.
    # We catch struct.error (or StructError) per-entry so that limited
    # environments like Brython's VFS _struct.py (which may not implement
    # every format code, e.g. '?') don't abort the whole table build.
    _StructError = getattr(_struct, 'error',
                   getattr(_struct, 'StructError', Exception))

    scalars = [
        ('b', 'b'), ('B', 'B'),
        ('h', 'h'), ('H', 'H'),
        ('i', 'i'), ('I', 'I'),
        ('q', 'q'), ('Q', 'Q'),
        ('f', 'f'), ('d', 'd'),
        ('?', '?'), ('c', 'c'),
    ]
    for code, fmt in scalars:
        try:
            sz = _struct.calcsize('=' + fmt)
        except (_StructError, Exception):
            # '?' (bool) unsupported by Brython _struct → treat as 1-byte int
            _FALLBACKS = {'?': ('B', 1, 1), 'c': ('B', 1, 1)}
            if code in _FALLBACKS:
                _TYPE_TABLE[code] = _FALLBACKS[code]
            continue
        try:
            native = _struct.calcsize(fmt)
        except (_StructError, Exception):
            native = sz
        aln = native if native <= sz else sz
        _TYPE_TABLE[code] = (fmt, sz, aln)

    # long / ulong: use native size (struct without '=' gives C ABI sizes)
    # e.g. calcsize('l')=8 on 64-bit Linux, =4 on Win32 or 32-bit
    for code, fmt in (('l', 'l'), ('L', 'L')):
        try:
            sz = _struct.calcsize(fmt)
        except (_StructError, Exception):
            sz = 4                              # safe fallback: 32-bit long
        aln      = min(sz, _PTR_SZ)
        pack_fmt = ('q' if sz == 8 else 'i') if code == 'l' else ('Q' if sz == 8 else 'I')
        _TYPE_TABLE[code] = (pack_fmt, sz, aln)

    # Pointer-width types — Brython is always 32-bit-ish JS; _PTR_SZ may be 4
    for code in ('P', 'z', 'Z', 'O'):
        _TYPE_TABLE[code] = (_PTR_FMT, _PTR_SZ, _PTR_SZ)

    # wchar_t: 2 bytes on Windows, 4 on Linux/Mac; Brython → 2
    wsz = 2 if _sys.platform in ('win32', 'brython') else 4
    _TYPE_TABLE['u'] = ('H' if wsz == 2 else 'I', wsz, wsz)

    # long double: may be 10 or 16 bytes; degrade to double for packing
    try:
        sz = _struct.calcsize('g')
        _TYPE_TABLE['g'] = ('d', sz, sz)
    except (_StructError, Exception):
        _TYPE_TABLE['g'] = ('d', 8, 8)


_build_type_table()


def _type_info(code: str) -> tuple[str, int, int]:
    if code not in _TYPE_TABLE:
        raise TypeError(f'unknown _type_ code {code!r}')
    return _TYPE_TABLE[code]


# ---------------------------------------------------------------------------
# Raw memory: array.array('B') as backing store
#
# array.array('B', ...).buffer_info() → (real_C_pointer_as_int, element_count)
# struct.pack_into / struct.unpack_from accept array objects directly.
# ---------------------------------------------------------------------------

class _RawBuf:
    """
    Fixed-size byte buffer backed by array.array('B').
    .address  → the real OS virtual address of the first byte.
    .size     → buffer size in bytes.
    Supports read(offset, n), write(offset, bytes), pack_into, unpack_from.
    """
    __slots__ = ('_arr', '_size', 'address')

    def __init__(self, size: int, init: bytes | None = None):
        if size == 0:
            size = 1                      # avoid zero-length edge cases
        self._arr  = _array.array('B', bytes(size))
        self._size = size
        self.address, _ = self._arr.buffer_info()
        if init:
            n = min(len(init), size)
            self._arr[:n] = _array.array('B', init[:n])

    # ------------------------------------------------------------------
    def read(self, offset: int = 0, size: int | None = None) -> bytes:
        end = self._size if size is None else offset + size
        return bytes(self._arr[offset:end])

    def write(self, data: _typing.Union[bytes, bytearray, _array.array], offset: int = 0):
        if isinstance(data, _array.array):
            data = bytes(data)
        self._arr[offset:offset + len(data)] = _array.array('B', data)

    def pack_into(self, fmt: str, offset: int, *values):
        _struct.pack_into(_ENDIAN + fmt, self._arr, offset, *values)

    def unpack_from(self, fmt: str, offset: int):
        return _struct.unpack_from(_ENDIAN + fmt, self._arr, offset)

    # ------------------------------------------------------------------
    def resize(self, new_size: int):
        """Grow the buffer, preserving existing content. Address may change."""
        if new_size <= self._size:
            return
        new_arr = _array.array('B', bytes(new_size))
        new_arr[:self._size] = self._arr
        self._arr  = new_arr
        self._size = new_size
        self.address, _ = self._arr.buffer_info()
        # Update global registry
        _addr_to_buf[self.address] = self

    def __len__(self) -> int:
        return self._size


# Global address registry: real_address (int) → _RawBuf
_addr_to_buf: dict[int, _RawBuf] = {}


def _alloc(size: int, init: bytes | None = None) -> _RawBuf:
    buf = _RawBuf(size, init)
    _addr_to_buf[buf.address] = buf
    return buf


def _buf_from_addr(addr: int) -> _RawBuf | None:
    return _addr_to_buf.get(addr)


def _buf_containing_addr(addr: int) -> _RawBuf | None:
    """Return the _RawBuf that *contains* addr (not just starts at it)."""
    # Fast path: exact match
    b = _addr_to_buf.get(addr)
    if b is not None:
        return b
    # Slow path: linear scan for the buffer spanning [base, base+size)
    for base, buf in _addr_to_buf.items():
        if base <= addr < base + len(buf):
            return buf
    return None


# ---------------------------------------------------------------------------
# sizeof / alignment  (accept type or instance)
# ---------------------------------------------------------------------------

def sizeof(obj) -> int:
    tp = obj if isinstance(obj, type) else type(obj)
    if hasattr(tp, '_size_'):
        return tp._size_
    raise TypeError(f'this type has no size: {tp!r}')


def alignment(obj) -> int:
    tp = obj if isinstance(obj, type) else type(obj)
    if hasattr(tp, '_align_'):
        return tp._align_
    raise TypeError(f'this type has no alignment: {tp!r}')


# ---------------------------------------------------------------------------
# addressof / byref / buffer_info / resize / cast / pointer
# ---------------------------------------------------------------------------

def addressof(obj) -> int:
    buf = getattr(obj, '_buffer_', None)
    if buf is None:
        raise TypeError('invalid type')
    return buf.address + getattr(obj, '_offset_', 0)


class _CArgObject:
    """Returned by byref() — a lightweight reference to a ctypes instance."""
    __slots__ = ('_obj', '_offset')

    def __init__(self, obj, offset: int = 0):
        self._obj    = obj
        self._offset = offset


def byref(obj, offset: int = 0) -> _CArgObject:
    if not isinstance(obj, _CData):
        raise TypeError(
            'byref() argument must be a ctypes instance, not %r'
            % type(obj).__name__)
    return _CArgObject(obj, offset)


def buffer_info(obj) -> tuple[int, int]:
    """Return (address, size_in_bytes)."""
    return addressof(obj), sizeof(obj)


def resize(obj, size: int):
    if not isinstance(obj, _CData):
        raise TypeError('expected ctypes instance')
    min_size = sizeof(obj)
    if size < min_size:
        raise ValueError('minimum size is %d' % min_size)
    old_addr = obj._buffer_.address
    obj._buffer_.resize(size)
    # If address changed, re-register
    new_addr = obj._buffer_.address
    if old_addr != new_addr:
        _addr_to_buf.pop(old_addr, None)
        _addr_to_buf[new_addr] = obj._buffer_


def pointer(obj) -> '_Pointer':
    """Create a typed pointer to obj."""
    pt   = POINTER(type(obj))
    inst = pt.__new__(pt)
    inst._buffer_       = _alloc(pt._size_)
    inst._offset_       = 0
    inst._b_needsfree_  = True
    inst._objects       = {id(obj): obj}
    inst._buffer_.pack_into(_PTR_FMT, 0, addressof(obj))
    return inst


def cast(obj, typ):
    """Reinterpret obj's address as a different type."""
    if isinstance(obj, _CArgObject):
        addr = addressof(obj._obj) + obj._offset
    elif isinstance(obj, _CData):
        addr = addressof(obj)
    elif isinstance(obj, int):
        addr = obj
    else:
        raise TypeError('cast() argument must be a pointer or integer')

    if isinstance(typ, type) and issubclass(typ, _Pointer):
        inst = typ.__new__(typ)
        inst._buffer_      = _alloc(typ._size_)
        inst._offset_      = 0
        inst._b_needsfree_ = True
        inst._objects      = {}
        inst._buffer_.pack_into(_PTR_FMT, 0, addr)
        return inst

    return typ.from_address(addr)


# ---------------------------------------------------------------------------
# Base _CData
# ---------------------------------------------------------------------------

class _CData:
    """Abstract base for all ctypes instances."""
    _b_needsfree_ = False
    _b_base_       = None
    _objects       = None
    _buffer_: _RawBuf | None = None
    _offset_: int = 0

    def __init_subclass__(cls, **kw):
        super().__init_subclass__(**kw)


# ---------------------------------------------------------------------------
# PyCSimpleType  →  _SimpleCData
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
        if not isinstance(source, (bytearray, memoryview)):
            raise TypeError('a bytes-like object is required')
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
        """Convert a Python value to a ctypes instance suitable for a function call."""
        if isinstance(value, cls):
            return value
        # For pointer types (z, Z, P) allow None → NULL
        if hasattr(cls, '_type_') and cls._type_ in ('z', 'Z', 'P') and value is None:
            return cls()
        try:
            return cls(value)
        except (TypeError, ValueError) as e:
            raise TypeError(
                f"wrong type for argument: expected {cls.__name__}, got {type(value).__name__}"
            ) from e


class _ExternalBuf:
    """
    Thin view over external (unowned) memory at a known address.
    Reads and writes use the global _addr_to_buf registry; if the address
    was allocated by us we delegate there, otherwise operations are no-ops
    (we cannot dereference arbitrary OS addresses without ctypes/cffi).
    """
    __slots__ = ('address', '_size')

    def __init__(self, address: int, size: int):
        self.address = address
        self._size   = size

    def _own(self) -> _RawBuf | None:
        return _addr_to_buf.get(self.address)

    def read(self, offset: int = 0, size: int | None = None) -> bytes:
        own = self._own()
        if own:
            return own.read(offset + (self.address - own.address), size)
        return bytes(size or self._size)            # zeros for unknown memory

    def write(self, data, offset: int = 0):
        own = self._own()
        if own:
            own.write(data, offset + (self.address - own.address))

    def pack_into(self, fmt: str, offset: int, *values):
        own = self._own()
        if own:
            own.pack_into(fmt, offset + (self.address - own.address), *values)

    def unpack_from(self, fmt: str, offset: int):
        own = self._own()
        if own:
            return own.unpack_from(fmt, offset + (self.address - own.address))
        sz = _struct.calcsize(fmt)
        return _struct.unpack(_ENDIAN + fmt, bytes(sz))

    def resize(self, new_size: int):
        self._size = new_size

    def __len__(self) -> int:
        return self._size


class _SimpleCData(_CData, metaclass=PyCSimpleType):
    """Base class for all simple ctypes scalar types."""

    def __init__(self, value=None):
        self._buffer_      = _alloc(self._size_)
        self._offset_      = 0
        self._b_needsfree_ = True
        if value is not None:
            self.value = value

    # ------------------------------------------------------------------
    @property
    def value(self):
        return self._read_value()

    @value.setter
    def value(self, v):
        self._write_value(v)

    # ------------------------------------------------------------------
    def _read_value(self):
        code = self._type_
        fmt  = self._fmt_
        off  = self._offset_

        if code == 'c':
            return self._buffer_.read(off, 1)

        if code == 'u':
            wsz = self._size_
            raw = self._buffer_.read(off, wsz)
            enc = 'utf-16-le' if wsz == 2 else 'utf-32-le'
            return raw.decode(enc)

        if code == 'z':
            addr = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            if addr == 0:
                return None
            return _read_cstring(addr)

        if code == 'Z':
            addr = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            if addr == 0:
                return None
            return _read_cwstring(addr)

        if code == 'O':
            # py_object: stored as an index into the _py_object_store dict
            idx = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            return _py_object_store.get(idx)

        return self._buffer_.unpack_from(fmt, off)[0]

    # ------------------------------------------------------------------
    def _write_value(self, v):
        code = self._type_
        fmt  = self._fmt_
        off  = self._offset_

        if code == 'c':
            if isinstance(v, int):
                v = bytes([v & 0xFF])
            if not isinstance(v, (bytes, bytearray)) or len(v) != 1:
                raise TypeError(
                    'one character bytes, bytearray or integer expected')
            self._buffer_.write(v, off)
            return

        if code == 'u':
            if not isinstance(v, str) or len(v) != 1:
                raise TypeError('a unicode character is required')
            enc = 'utf-16-le' if self._size_ == 2 else 'utf-32-le'
            self._buffer_.write(v.encode(enc), off)
            return

        if code in ('z', 'Z', 'P'):
            if v is None:
                self._buffer_.pack_into(_PTR_FMT, off, 0)
            elif isinstance(v, int):
                self._buffer_.pack_into(_PTR_FMT, off, v)
            elif isinstance(v, (bytes, bytearray)):
                s = bytes(v) + b'\x00'
                tmp = _alloc(len(s), s)
                self._buffer_.pack_into(_PTR_FMT, off, tmp.address)
                if self._objects is None:
                    self._objects = {}
                self._objects[tmp.address] = tmp
            elif isinstance(v, str):
                s = v.encode() + b'\x00'
                tmp = _alloc(len(s), s)
                self._buffer_.pack_into(_PTR_FMT, off, tmp.address)
                if self._objects is None:
                    self._objects = {}
                self._objects[tmp.address] = tmp
            return

        if code == 'O':
            idx = id(v)
            _py_object_store[idx] = v
            self._buffer_.pack_into(_PTR_FMT, off, idx)
            return

        if code == '?':
            v = int(bool(v))

        try:
            self._buffer_.pack_into(fmt, off, v)
        except _struct.error as e:
            raise OverflowError(str(e)) from e

    # ------------------------------------------------------------------
    def __repr__(self):
        try:
            return f'{type(self).__name__}({self.value!r})'
        except Exception:
            return f'{type(self).__name__}(<error>)'

    def __bool__(self):
        try:
            return bool(self.value)
        except Exception:
            return False

    def __eq__(self, other):
        if isinstance(other, _SimpleCData):
            return self.value == other.value
        return self.value == other

    def __hash__(self):
        return hash(self.value)


# py_object backing store (maps id(obj) → obj, preventing GC)
_py_object_store: dict[int, object] = {}

# Null-terminated string readers (work only for our own heap)
def _read_cstring(addr: int) -> bytes:
    buf = _buf_from_addr(addr)
    if buf is None:
        return b''
    raw = buf.read(addr - buf.address)
    nul = raw.find(b'\x00')
    return raw[:nul] if nul >= 0 else raw


def _read_cwstring(addr: int) -> str:
    buf = _buf_from_addr(addr)
    if buf is None:
        return ''
    wsz  = 4 if _sys.platform != 'win32' else 2
    raw  = buf.read(addr - buf.address)
    enc  = 'utf-32-le' if wsz == 4 else 'utf-16-le'
    # find double-null (or quad-null for utf-32)
    for i in range(0, len(raw) - wsz + 1, wsz):
        if all(b == 0 for b in raw[i:i + wsz]):
            return raw[:i].decode(enc)
    return raw.decode(enc, errors='replace')


# ---------------------------------------------------------------------------
# Array type factory
# ---------------------------------------------------------------------------

def _make_array_type(item_type, count: int):
    key = (item_type, count)
    cached = _pointer_type_cache.get(key)
    if cached is not None:
        return cached

    sz   = sizeof(item_type) * count
    aln  = alignment(item_type)
    name = f'{item_type.__name__}_Array_{count}'

    ns = {
        '_type_':   item_type,
        '_length_': count,
        '_size_':   sz,
        '_align_':  aln,
    }
    cls = PyCArrayType(name, (Array,), ns)
    _pointer_type_cache[key] = cls
    return cls


class PyCArrayType(type):

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


class Array(_CData, metaclass=PyCArrayType):
    """Base class for ctypes array types."""

    def __init__(self, *args):
        self._buffer_      = _alloc(self._size_)
        self._offset_      = 0
        self._b_needsfree_ = True
        for i, v in enumerate(args):
            if i >= self._length_:
                break
            self._write_item(i, v)

    # ------------------------------------------------------------------
    def _item_offset(self, index: int) -> int:
        return self._offset_ + index * sizeof(self._type_)

    def _write_item(self, index: int, value):
        item_sz = sizeof(self._type_)
        off     = self._item_offset(index)
        if issubclass(self._type_, _SimpleCData):
            tmp = self._type_(value)
            self._buffer_.write(tmp._buffer_.read(0, item_sz), off)
        elif isinstance(value, self._type_):
            self._buffer_.write(value._buffer_.read(0, item_sz), off)
        elif isinstance(value, _CData):
            self._buffer_.write(value._buffer_.read(0, item_sz), off)
        else:
            # Try coercion
            tmp = self._type_(value)
            self._buffer_.write(tmp._buffer_.read(0, item_sz), off)

    def _read_item(self, index: int):
        off  = self._item_offset(index)
        addr = self._buffer_.address + off
        inst = self._type_.__new__(self._type_)
        inst._buffer_      = self._buffer_
        inst._offset_      = off
        inst._b_needsfree_ = False
        return inst

    # ------------------------------------------------------------------
    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(self._length_))]
        if index < 0:
            index += self._length_
        if not 0 <= index < self._length_:
            raise IndexError('array index out of range')
        inst = self._read_item(index)
        if issubclass(self._type_, _SimpleCData):
            return inst.value
        return inst

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            for i, v in zip(range(*index.indices(self._length_)), value):
                self._write_item(i, v)
            return
        if index < 0:
            index += self._length_
        if not 0 <= index < self._length_:
            raise IndexError('array index out of range')
        self._write_item(index, value)

    def __len__(self) -> int:
        return self._length_

    def __iter__(self):
        for i in range(self._length_):
            yield self[i]

    # ------------------------------------------------------------------
    @property
    def value(self):
        if issubclass(self._type_, _SimpleCData):
            code = self._type_._type_
            if code == 'c':
                raw = self._buffer_.read(self._offset_, self._size_)
                nul = raw.find(b'\x00')
                return raw[:nul] if nul >= 0 else raw
            if code == 'u':
                wsz = self._type_._size_
                raw = self._buffer_.read(self._offset_, self._size_)
                enc = 'utf-16-le' if wsz == 2 else 'utf-32-le'
                return raw.decode(enc, errors='replace').rstrip('\x00')
        raise AttributeError(f'{type(self).__name__} has no .value')

    @value.setter
    def value(self, v):
        if issubclass(self._type_, _SimpleCData):
            code = self._type_._type_
            if code == 'c':
                if isinstance(v, str):
                    v = v.encode()
                pad = self._size_ - len(v) - 1
                raw = (v[:self._size_ - 1] + b'\x00' +
                       b'\x00' * max(pad, 0))
                self._buffer_.write(raw[:self._size_], self._offset_)
                return
            if code == 'u':
                wsz = self._type_._size_
                enc = 'utf-16-le' if wsz == 2 else 'utf-32-le'
                raw = v.encode(enc)
                cap = self._size_ - wsz
                raw = raw[:cap] + b'\x00' * wsz
                self._buffer_.write(raw.ljust(self._size_, b'\x00'),
                                    self._offset_)
                return
        raise AttributeError(f'cannot set .value on {type(self).__name__}')

    def __repr__(self):
        return f'<{type(self).__name__} object at 0x{addressof(self):016x}>'

    # ------------------------------------------------------------------
    @property
    def raw(self):
        """Raw bytes of this array (c_char arrays only in real ctypes,
        but we allow it on any array for memmove / memset compat)."""
        return self._buffer_.read(self._offset_, self._size_)

    @raw.setter
    def raw(self, data: bytes):
        if len(data) > self._size_:
            raise ValueError('raw buffer is too large')
        self._buffer_.write(data.ljust(self._size_, b'\x00'), self._offset_)


# ---------------------------------------------------------------------------
# Structure / Union field layout
# ---------------------------------------------------------------------------

class _FieldDescriptor:
    """Data descriptor for a single field in a Structure or Union."""

    __slots__ = ('name', 'ctype', 'offset', 'size',
                 'bit_offset', 'bit_size')

    def __init__(self, name, ctype, offset, size,
                 bit_offset=0, bit_size=0):
        self.name       = name
        self.ctype      = ctype
        self.offset     = offset
        self.size       = size
        self.bit_offset = bit_offset
        self.bit_size   = bit_size

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        off = obj._offset_ + self.offset
        if self.bit_size:
            return self._get_bits(obj, off)
        if issubclass(self.ctype, _SimpleCData):
            inst = self.ctype.__new__(self.ctype)
            inst._buffer_      = obj._buffer_
            inst._offset_      = off
            inst._b_needsfree_ = False
            return inst.value
        # Composite / array / pointer: return a live view
        inst = self.ctype.__new__(self.ctype)
        inst._buffer_      = obj._buffer_
        inst._offset_      = off
        inst._b_needsfree_ = False
        return inst

    def __set__(self, obj, value):
        off = obj._offset_ + self.offset
        if self.bit_size:
            self._set_bits(obj, off, value)
            return
        if issubclass(self.ctype, _SimpleCData):
            tmp = self.ctype(value)
            obj._buffer_.write(tmp._buffer_.read(0, self.size), off)
        elif isinstance(value, self.ctype):
            obj._buffer_.write(value._buffer_.read(0, self.size), off)
        elif isinstance(value, _CData):
            obj._buffer_.write(value._buffer_.read(0, self.size), off)
        else:
            # Try coercion via the type's constructor
            tmp = self.ctype(value)
            obj._buffer_.write(tmp._buffer_.read(0, self.size), off)

    def _get_bits(self, obj, off):
        fmt  = self.ctype._fmt_
        word = obj._buffer_.unpack_from(fmt, off)[0]
        mask = (1 << self.bit_size) - 1
        return (word >> self.bit_offset) & mask

    def _set_bits(self, obj, off, value):
        fmt  = self.ctype._fmt_
        word = obj._buffer_.unpack_from(fmt, off)[0]
        mask = (1 << self.bit_size) - 1
        word = (word & ~(mask << self.bit_offset)) | \
               ((int(value) & mask) << self.bit_offset)
        obj._buffer_.pack_into(fmt, off, word)


def _layout_fields(fields, pack=None, is_union=False):
    """
    Compute field offsets, descriptors, total size, and max alignment.
    Returns (descriptors_dict, total_size, max_align).
    """
    offset    = 0
    max_align = 1
    descriptors = {}

    for spec in fields:
        fname    = spec[0]
        ctype    = spec[1]
        bit_size = spec[2] if len(spec) > 2 else 0

        sz  = sizeof(ctype)
        aln = alignment(ctype)
        if pack is not None:
            aln = min(aln, pack)
        max_align = max(max_align, aln)

        if is_union:
            field_off = 0
        else:
            if aln > 1 and offset % aln:
                offset += aln - (offset % aln)
            field_off = offset

        if bit_size:
            bit_off = 0 if is_union else (offset - field_off) * 8
            descriptors[fname] = _FieldDescriptor(
                fname, ctype, field_off, sz, bit_off, bit_size)
            if not is_union:
                offset = field_off + sz
        else:
            descriptors[fname] = _FieldDescriptor(
                fname, ctype, field_off, sz)
            if not is_union:
                offset = field_off + sz

    if is_union:
        total = max((sizeof(s[1]) for s in fields), default=0)
    else:
        total = offset

    # Pad total to struct alignment
    if max_align > 1 and total % max_align:
        total += max_align - (total % max_align)

    return descriptors, total, max_align


# ---------------------------------------------------------------------------
# PyCStructType  →  Structure
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# UnionType  →  Union
# ---------------------------------------------------------------------------

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
# PyCPointerType  →  _Pointer
# ---------------------------------------------------------------------------

def POINTER(ctype):
    """Return the pointer type for ctype, creating it if necessary."""
    key = ('ptr', ctype)
    cached = _pointer_type_cache.get(key)
    if cached is not None:
        return cached

    name = f'LP_{ctype.__name__}'
    ns   = {
        '_type_':  ctype,
        '_size_':  _PTR_SZ,
        '_align_': _PTR_SZ,
        '_fmt_':   _PTR_FMT,
    }
    cls = PyCPointerType(name, (_Pointer,), ns)
    _pointer_type_cache[key] = cls
    return cls


class PyCPointerType(type):

    def __mul__(cls, count: int):
        return _make_array_type(cls, count)

    def from_address(cls, address: int):
        buf  = _buf_containing_addr(address) or _ExternalBuf(address, cls._size_)
        inst = cls.__new__(cls)
        inst._buffer_      = buf
        inst._offset_      = address - buf.address
        inst._b_needsfree_ = False
        inst._objects      = {}
        return inst


class _Pointer(_CData, metaclass=PyCPointerType):
    """Base class for ctypes pointer types."""

    _size_  = _PTR_SZ
    _align_ = _PTR_SZ
    _fmt_   = _PTR_FMT

    def __init__(self, obj=None):
        self._buffer_      = _alloc(self._size_)
        self._offset_      = 0
        self._b_needsfree_ = True
        self._objects      = {}
        if obj is not None:
            addr = addressof(obj)
            self._buffer_.pack_into(self._fmt_, self._offset_, addr)
            self._objects[id(obj)] = obj

    # ------------------------------------------------------------------
    def _get_ptr_val(self) -> int:
        return self._buffer_.unpack_from(self._fmt_, self._offset_)[0]

    def _set_ptr_val(self, addr: int):
        self._buffer_.pack_into(self._fmt_, self._offset_, addr)

    # ------------------------------------------------------------------
    @property
    def contents(self):
        addr = self._get_ptr_val()
        if addr == 0:
            raise ValueError('NULL pointer access')
        return self._type_.from_address(addr)

    @contents.setter
    def contents(self, obj):
        addr = addressof(obj)
        self._set_ptr_val(addr)
        self._objects[id(obj)] = obj

    def __getitem__(self, index: int):
        addr = self._get_ptr_val()
        if addr == 0:
            raise ValueError('NULL pointer dereference')
        item_sz = sizeof(self._type_)
        target  = addr + index * item_sz
        inst    = self._type_.from_address(target)
        if issubclass(self._type_, _SimpleCData):
            return inst.value
        return inst

    def __setitem__(self, index: int, value):
        addr = self._get_ptr_val()
        if addr == 0:
            raise ValueError('NULL pointer dereference')
        item_sz = sizeof(self._type_)
        target  = addr + index * item_sz
        buf     = _buf_from_addr(addr)
        if buf is None:
            raise ValueError('pointer points outside managed heap')
        if issubclass(self._type_, _SimpleCData):
            tmp  = self._type_(value)
            data = tmp._buffer_.read(0, item_sz)
            off  = target - buf.address
            buf.write(data, off)
        else:
            raise TypeError(
                'pointer item assignment not supported for this type')

    @property
    def value(self) -> int:
        return self._get_ptr_val()

    def __bool__(self) -> bool:
        return self._get_ptr_val() != 0

    def __repr__(self):
        return (f'<{type(self).__name__} at 0x{addressof(self):016x}'
                f' → 0x{self._get_ptr_val():016x}>')


# ---------------------------------------------------------------------------
# Virtual library / FFI layer
#
# Without ctypes/cffi we cannot call native C functions.
# Applications register Python callables as "native functions" via
# register_library() / register_symbol(), and the normal CFuncPtr
# machinery invokes them.
# ---------------------------------------------------------------------------

# _vlibs: { handle_int: { symbol_name: (callable, address_int) } }
# Handle 0 is reserved for dlopen(None) — the current process / Python C API.
# Symbols looked up on it get synthetic addresses; calling them raises unless
# a real callable has been registered via register_library().
_PROC_HANDLE = 0
_vlibs:          dict[int, dict[str, tuple]] = {_PROC_HANDLE: {}}
_vlib_names:     dict[int, str]              = {_PROC_HANDLE: '<process>'}
_vcallables:     dict[int, callable]         = {}  # addr → Python callable
_vlib_counter    = 0
_vaddr_counter   = 0x7FFF_0000_0000_0000     # high synthetic address space


def _next_vaddr() -> int:
    global _vaddr_counter
    _vaddr_counter += 8
    return _vaddr_counter


def register_library(name: str, symbols: dict[str, callable]) -> int:
    """
    Register a Python-callable library under a name.

    Usage::

        import math
        handle = register_library('libm.so', {
            'cos': math.cos,
            'sin': math.sin,
        })

    Returns an opaque integer handle for use with dlsym().
    """
    global _vlib_counter
    _vlib_counter += 1
    h    = _vlib_counter
    syms = {}
    for sym_name, callable_ in symbols.items():
        addr = _next_vaddr()
        syms[sym_name] = (callable_, addr)
        _vcallables[addr] = callable_
    _vlibs[h]      = syms
    _vlib_names[h] = name
    return h


def _vlib_sym(lib_handle: int, name: str) -> int:
    """Return synthetic address for a registered symbol.
    For the process handle (0 / dlopen(None)), unknown symbols get a
    freshly-minted virtual address so lookups never raise — calls will
    raise at call-time if no real callable is registered.
    """
    syms = _vlibs.get(lib_handle)
    if syms is None:
        raise OSError(f'invalid library handle {lib_handle}')
    entry = syms.get(name)
    if entry is None:
        if lib_handle == _PROC_HANDLE:
            # Auto-allocate a stub address; calling it will raise OSError.
            addr = _next_vaddr()
            syms[name] = (None, addr)
            return addr
        raise OSError(f'symbol {name!r} not found in {_vlib_names[lib_handle]!r}')
    return entry[1]


def dlopen(name: str | None, mode: int = RTLD_LOCAL) -> int:
    """
    Open a library by name.

    ``None`` means "the current process" (same as POSIX dlopen(NULL, ...)).
    Returns the reserved process handle so callers like ``PyDLL(None)`` succeed.

    For any other name, if it matches a previously registered library name,
    returns its handle.  Otherwise raises OSError — use register_library() to
    provide implementations.
    """
    if name is None:
        return _PROC_HANDLE
    for h, n in _vlib_names.items():
        if n == name or _os.path.basename(n) == _os.path.basename(name):
            return h
    raise OSError(
        f'dlopen({name!r}): no registered implementation — '
        f'use register_library() to provide one')


def dlclose(handle: int) -> None:
    _vlibs.pop(handle, None)
    _vlib_names.pop(handle, None)


def dlsym(handle: int, name: str) -> int:
    return _vlib_sym(handle, name)


# ---------------------------------------------------------------------------
# PyCFuncPtrType  →  CFuncPtr
# ---------------------------------------------------------------------------

class PyCFuncPtrType(type):

    def __mul__(cls, count: int):
        return _make_array_type(cls, count)


class CFuncPtr(_CData, metaclass=PyCFuncPtrType):
    """Base class for ctypes function pointer types."""

    _size_  = _PTR_SZ
    _align_ = _PTR_SZ
    restype  = None
    argtypes = None
    errcheck = None
    _flags_  = FUNCFLAG_CDECL

    def __init__(self, *args):
        self._buffer_      = _alloc(self._size_)
        self._offset_      = 0
        self._b_needsfree_ = True
        self._callable_    = None
        self._func_addr_   = 0

        if not args:
            return

        first = args[0]

        if callable(first) and not isinstance(first, int):
            # Python callback
            self._callable_ = first
            addr = _next_vaddr()
            _vcallables[addr] = first
            self._func_addr_  = addr
            self._buffer_.pack_into(_PTR_FMT, 0, addr)

        elif isinstance(first, int):
            self._func_addr_ = first
            self._buffer_.pack_into(_PTR_FMT, 0, first)

        elif isinstance(first, tuple) and len(first) == 2:
            name_or_ord, lib = first
            if isinstance(name_or_ord, str):
                lib_handle = lib if isinstance(lib, int) else getattr(lib, '_handle', id(lib))
                addr = _vlib_sym(lib_handle, name_or_ord)
            else:
                addr = int(name_or_ord)
            self._func_addr_ = addr
            self._buffer_.pack_into(_PTR_FMT, 0, addr)

    # ------------------------------------------------------------------
    def __call__(self, *args):
        addr = self._buffer_.unpack_from(_PTR_FMT, self._offset_)[0]
        if addr == 0:
            raise OSError('call through NULL function pointer')

        cb = self._callable_ or _vcallables.get(addr)
        if cb is None:
            raise OSError(
                f'no Python callable registered for address 0x{addr:x} — '
                f'use register_library() / CFuncPtr(callback) to provide one')

        # Unwrap ctypes arguments to Python values
        py_args = []
        for i, arg in enumerate(args):
            if isinstance(arg, _CArgObject):
                py_args.append(arg._obj)
            elif isinstance(arg, _SimpleCData):
                py_args.append(arg.value)
            else:
                py_args.append(arg)

        result = cb(*py_args)

        # Optionally wrap result in restype
        if self.errcheck is not None:
            result = self.errcheck(result, self, args)

        if self.restype is not None and not isinstance(result, self.restype):
            try:
                result = self.restype(result).value
            except Exception:
                pass

        return result

    def __repr__(self):
        return f'<{type(self).__name__} at 0x{addressof(self):016x}>'


# ---------------------------------------------------------------------------
# Low-level call functions
# ---------------------------------------------------------------------------

def call_function(func_addr: int, arguments: tuple):
    """Call the function at func_addr (stdcall/cdecl)."""
    cb = _vcallables.get(func_addr)
    if cb is None:
        raise OSError(
            f'call_function: no callable at 0x{func_addr:x} — '
            f'register it via register_library()')
    return cb(*arguments)


def call_cdeclfunction(func_addr: int, arguments: tuple):
    """Call a cdecl function at func_addr."""
    return call_function(func_addr, arguments)


# ---------------------------------------------------------------------------
# errno (virtual)
# ---------------------------------------------------------------------------
_errno_value = 0

def get_errno() -> int:
    return _errno_value

def set_errno(value: int) -> None:
    global _errno_value
    _errno_value = int(value)


# ---------------------------------------------------------------------------
# PyObj_FromPtr / Py_INCREF / Py_DECREF
# ---------------------------------------------------------------------------

def PyObj_FromPtr(addr: int):
    """Recover a Python object stored in the py_object backing store."""
    return _py_object_store.get(addr)


def Py_INCREF(obj):
    pass   # handled by Python GC


def Py_DECREF(obj):
    pass   # handled by Python GC


# ---------------------------------------------------------------------------
# Special addresses — each gets a real Python callable wired into _vcallables
# so that CFuncPtr.__call__ can dispatch them correctly.
# ctypes/__init__.py calls: _cast(obj, obj, typ)  (obj passed twice)
# ---------------------------------------------------------------------------

def _impl_cast(obj, obj2, typ):
    """Backend for ctypes.cast() routed via _cast_addr."""
    return cast(obj2 if obj2 is not None else obj, typ)

def _impl_memmove(dst, src, count):
    """Backend for ctypes.memmove()."""
    # dst / src may be ints (addresses) or ctypes instances
    def _addr(x):
        if isinstance(x, int):     return x
        if isinstance(x, _CData):  return addressof(x)
        if isinstance(x, _CArgObject): return addressof(x._obj) + x._offset
        raise TypeError(f'memmove: bad argument {type(x)}')

    dst_addr = _addr(dst)
    src_addr = _addr(src)
    n        = int(count)

    src_buf = _buf_from_addr(src_addr)
    dst_buf = _buf_from_addr(dst_addr)

    if src_buf is not None:
        data = src_buf.read(src_addr - src_buf.address, n)
    else:
        data = bytes(n)             # unknown src → zeros

    if dst_buf is not None:
        dst_buf.write(data, dst_addr - dst_buf.address)
    # if dst is unknown external memory we silently skip (can't write there)

    return dst_addr

def _impl_memset(dst, c, count):
    """Backend for ctypes.memset()."""
    def _addr(x):
        if isinstance(x, int):    return x
        if isinstance(x, _CData): return addressof(x)
        raise TypeError(f'memset: bad argument {type(x)}')

    dst_addr = _addr(dst)
    n        = int(count)
    val      = int(c) & 0xFF

    dst_buf = _buf_from_addr(dst_addr)
    if dst_buf is not None:
        dst_buf.write(bytes([val] * n), dst_addr - dst_buf.address)

    return dst_addr

def _impl_string_at(addr, size=-1):
    """Backend for ctypes.string_at()."""
    if isinstance(addr, _CData):
        addr = addressof(addr)
    addr = int(addr)
    buf  = _buf_from_addr(addr)
    if buf is None:
        return b''
    offset = addr - buf.address
    raw    = buf.read(offset)
    if size < 0:
        # null-terminated
        nul = raw.find(b'\x00')
        return raw[:nul] if nul >= 0 else raw
    return raw[:size]

def _impl_wstring_at(addr, size=-1):
    """Backend for ctypes.wstring_at()."""
    if isinstance(addr, _CData):
        addr = addressof(addr)
    addr  = int(addr)
    wsz   = 2 if _sys.platform == 'win32' else 4
    enc   = 'utf-16-le' if wsz == 2 else 'utf-32-le'
    buf   = _buf_from_addr(addr)
    if buf is None:
        return ''
    offset = addr - buf.address
    raw    = buf.read(offset)
    if size < 0:
        # find null wchar
        i = 0
        while i + wsz <= len(raw):
            if raw[i:i+wsz] == b'\x00' * wsz:
                break
            i += wsz
        raw = raw[:i]
    else:
        raw = raw[:size * wsz]
    return raw.decode(enc, errors='replace')

_cast_addr       = _next_vaddr()
_memmove_addr    = _next_vaddr()
_memset_addr     = _next_vaddr()
_string_at_addr  = _next_vaddr()
_wstring_at_addr = _next_vaddr()

# Register implementations so CFuncPtr.__call__ can find them
_vcallables[_cast_addr]       = _impl_cast
_vcallables[_memmove_addr]    = _impl_memmove
_vcallables[_memset_addr]     = _impl_memset
_vcallables[_string_at_addr]  = _impl_string_at
_vcallables[_wstring_at_addr] = _impl_wstring_at


# ---------------------------------------------------------------------------
# _unpickle
# ---------------------------------------------------------------------------

def _unpickle(typ, state: bytes):
    inst = typ.__new__(typ)
    inst._buffer_      = _alloc(sizeof(typ), state)
    inst._offset_      = 0
    inst._b_needsfree_ = True
    return inst


# ---------------------------------------------------------------------------
# Public API surface  (mirrors real _ctypes.__all__)
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
    # extensions
    'register_library',
]

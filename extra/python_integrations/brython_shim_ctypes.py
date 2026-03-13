from __future__ import annotations
# Dylan Kenneth Eliot - Brython Optimized _ctypes shim

"""
This is the brython shim being tested on premise that it will work during testing & has worked for cpython similarly.

"""

"""
brython_shim_ctypes.py
======================
Fully-fleshed Brython shim for _ctypes.

Differences from _ctypes_pure.py (the CPython reference):
  • PTR_SZ is hard-wired to 4 (JS is a 32-bit address space for our virtual heap).
  • buffer_info() uses id() instead of array.buffer_info() (unavailable in Brython).
  • write() iterates byte-by-byte (Brython slice-assign on array may differ).
  • register_library() mirrors symbols onto window so JS code can call them.
  • JS bridge helpers (js_array_buffer, js_data_view, js_callback) let you pass
    ctypes buffers and Python callables into browser Web APIs.
  • All symbols declared in __all__ are actually defined below.

Dependencies: Brython standard library only + `from browser import window`.
"""

import array      as _array
import struct     as _struct
import sys        as _sys
import os         as _os
import typing     as _typing
import weakref
from collections import OrderedDict
from browser     import window

# ---------------------------------------------------------------------------
# Platform constants
# Derive pointer size from struct so _check_size(py_object, "P") passes
# regardless of whether the underlying JS engine reports 4 or 8 bytes.
# ---------------------------------------------------------------------------
_PTR_SZ  = _struct.calcsize('P')          # 8 on 64-bit JS, 4 on 32-bit
_PTR_FMT = 'Q' if _PTR_SZ == 8 else 'I'  # unsigned int wide enough for a pointer
_ENDIAN  = '<'

__version__            = '1.1.0'
CTYPES_MAX_ARGCOUNT    = 1024
FUNCFLAG_CDECL         = 0x1
FUNCFLAG_PYTHONAPI     = 0x4
FUNCFLAG_USE_ERRNO     = 0x8
FUNCFLAG_USE_LASTERROR = 0x10
RTLD_LOCAL             = 0
RTLD_GLOBAL            = 0x100
SIZEOF_TIME_T          = 4          # 32-bit JS context

_pointer_type_cache: dict = {}

# ---------------------------------------------------------------------------
# ArgumentError
# ---------------------------------------------------------------------------
class ArgumentError(Exception):
    pass

# ---------------------------------------------------------------------------
# Type-code table
# ---------------------------------------------------------------------------
_TYPE_TABLE: dict = {}

def _build_type_table():
    _StructError = getattr(_struct, 'error', getattr(_struct, 'StructError', Exception))

    scalars = [
        ('b','b'),('B','B'),
        ('h','h'),('H','H'),
        ('i','i'),('I','I'),
        ('q','q'),('Q','Q'),
        ('f','f'),('d','d'),
        ('?','?'),('c','c'),
    ]
    for code, fmt in scalars:
        try:
            sz  = _struct.calcsize('=' + fmt)
            aln = sz
            _TYPE_TABLE[code] = (fmt, sz, aln)
        except Exception:
            # Brython VFS _struct may not support every format code
            _TYPE_TABLE[code] = ('B', 1, 1)

    # long / ulong  (native width — treat as 4 bytes in 32-bit JS context)
    for code, signed in (('l', True), ('L', False)):
        fmt = 'i' if signed else 'I'
        _TYPE_TABLE[code] = (fmt, 4, 4)

    # Pointer-width types
    for code in ('P', 'z', 'Z', 'O'):
        _TYPE_TABLE[code] = (_PTR_FMT, _PTR_SZ, _PTR_SZ)

    # wchar_t → 2 bytes (JS / Windows convention)
    _TYPE_TABLE['u'] = ('H', 2, 2)

    # long double → degrade to double
    _TYPE_TABLE['g'] = ('d', 8, 8)

_build_type_table()

def _type_info(code: str):
    return _TYPE_TABLE.get(code, ('B', 1, 1))

# ---------------------------------------------------------------------------
# Raw buffer  — backed by array.array('B'), addressed via id()
# ---------------------------------------------------------------------------

class _RawBuf:
    """
    Fixed-size byte buffer.  .address is id(self._arr), giving a
    unique integer that plays the role of a virtual memory address
    inside our managed heap.
    """
    __slots__ = ('_arr', '_size', 'address')

    def __init__(self, size: int, init=None):
        if size == 0:
            size = 1
        self._arr  = _array.array('B', [0] * size)
        self._size = size
        self.address = id(self._arr)     # Brython: id() as virtual address

        if init:
            n = min(len(init), size)
            for i in range(n):
                self._arr[i] = init[i] if isinstance(init[i], int) else init[i]

    # ------------------------------------------------------------------
    def read(self, offset: int = 0, size=None) -> bytes:
        end = self._size if size is None else min(offset + size, self._size)
        return bytes(self._arr[offset:end])

    def write(self, data, offset: int = 0):
        if isinstance(data, _array.array):
            data = bytes(data)
        for i, b in enumerate(data):
            pos = offset + i
            if pos < self._size:
                self._arr[pos] = b if isinstance(b, int) else ord(b)

    def pack_into(self, fmt: str, offset: int, *values):
        full_fmt = _ENDIAN + fmt
        # Pack to a plain bytes object first — avoids Brython VFS _struct
        # bug where pack_into does buf[o:o+n]=data on an array.array and
        # raises IndexError even after the array has been extended.
        data = _struct.pack(full_fmt, *values)
        needed = offset + len(data)
        if needed > self._size:
            self._arr.extend([0] * (needed - self._size))
            self._size = len(self._arr)
        for i, b in enumerate(data):
            self._arr[offset + i] = b if isinstance(b, int) else ord(b)

    def unpack_from(self, fmt: str, offset: int):
        full_fmt = _ENDIAN + fmt
        sz  = _struct.calcsize(full_fmt)
        end = min(offset + sz, self._size)
        raw = bytes(self._arr[offset:end]).ljust(sz, b'\x00')
        return _struct.unpack(full_fmt, raw)

    def resize(self, new_size: int):
        if new_size <= self._size:
            return
        self._arr.extend([0] * (new_size - self._size))
        self._size = new_size
        # address (id) does NOT change after extend — same object

    def __len__(self) -> int:
        return self._size


# Global virtual-heap registry
_addr_to_buf: dict = {}

def _alloc(size: int, init=None) -> _RawBuf:
    buf = _RawBuf(size, init)
    _addr_to_buf[buf.address] = buf
    return buf

def _buf_from_addr(addr: int):
    return _addr_to_buf.get(addr)

def _buf_containing_addr(addr: int):
    b = _addr_to_buf.get(addr)
    if b is not None:
        return b
    for base, buf in _addr_to_buf.items():
        if base <= addr < base + len(buf):
            return buf
    return None

# ---------------------------------------------------------------------------
# _ExternalBuf — thin view over unowned (or JS-side) memory
# ---------------------------------------------------------------------------

class _ExternalBuf:
    """
    Represents an address that lives outside our managed heap.
    Reads return zeros; writes are silently dropped unless the address
    happens to be in _addr_to_buf.
    """
    __slots__ = ('address', '_size')

    def __init__(self, address: int, size: int):
        self.address = address
        self._size   = size

    def _own(self):
        return _addr_to_buf.get(self.address)

    def read(self, offset: int = 0, size=None) -> bytes:
        own = self._own()
        if own:
            return own.read(offset, size)
        return bytes(size or self._size)

    def write(self, data, offset: int = 0):
        own = self._own()
        if own:
            own.write(data, offset)

    def pack_into(self, fmt: str, offset: int, *values):
        own = self._own()
        if own:
            own.pack_into(fmt, offset, *values)

    def unpack_from(self, fmt: str, offset: int):
        own = self._own()
        if own:
            return own.unpack_from(fmt, offset)
        sz = _struct.calcsize(_ENDIAN + fmt)
        return _struct.unpack(_ENDIAN + fmt, bytes(sz))

    def resize(self, new_size: int):
        self._size = new_size

    def __len__(self) -> int:
        return self._size

# ---------------------------------------------------------------------------
# sizeof / alignment / addressof / byref / buffer_info / resize / pointer / cast
# ---------------------------------------------------------------------------

def sizeof(obj) -> int:
    tp = obj if isinstance(obj, type) else type(obj)
    return getattr(tp, '_size_', 0)

def alignment(obj) -> int:
    tp = obj if isinstance(obj, type) else type(obj)
    return getattr(tp, '_align_', 1)

def addressof(obj) -> int:
    buf = getattr(obj, '_buffer_', None)
    if buf is None:
        raise TypeError('invalid type for addressof')
    return buf.address + getattr(obj, '_offset_', 0)

class _CArgObject:
    __slots__ = ('_obj', '_offset')
    def __init__(self, obj, offset: int = 0):
        self._obj    = obj
        self._offset = offset

def byref(obj, offset: int = 0) -> _CArgObject:
    if not isinstance(obj, _CData):
        raise TypeError('byref() argument must be a ctypes instance')
    return _CArgObject(obj, offset)

def buffer_info(obj):
    return addressof(obj), sizeof(obj)

def resize(obj, size: int):
    if not isinstance(obj, _CData):
        raise TypeError('expected ctypes instance')
    min_size = sizeof(obj)
    if size < min_size:
        raise ValueError('minimum size is %d' % min_size)
    obj._buffer_.resize(size)

def pointer(obj):
    pt   = POINTER(type(obj))
    inst = pt.__new__(pt)
    inst._buffer_      = _alloc(pt._size_)
    inst._offset_      = 0
    inst._b_needsfree_ = True
    inst._objects      = {id(obj): obj}
    inst._buffer_.pack_into(_PTR_FMT, 0, addressof(obj))
    return inst

def cast(obj, typ):
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
# py_object backing store + null-terminated string readers
# ---------------------------------------------------------------------------
_py_object_store: dict = {}

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
    wsz = 2          # always UTF-16 LE on Brython/JS
    enc = 'utf-16-le'
    raw = buf.read(addr - buf.address)
    for i in range(0, len(raw) - wsz + 1, wsz):
        if raw[i:i+wsz] == b'\x00\x00':
            return raw[:i].decode(enc)
    return raw.decode(enc, errors='replace')

# ---------------------------------------------------------------------------
# Base _CData
# ---------------------------------------------------------------------------

class _CData:
    _b_needsfree_ = False
    _b_base_       = None
    _objects       = None
    _buffer_       = None
    _offset_: int  = 0

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
        inst._buffer_      = buf
        inst._offset_      = address - buf.address
        inst._b_needsfree_ = False
        return inst

    def from_buffer(cls, source):
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
            raise TypeError(
                f'wrong type: expected {cls.__name__}, got {type(value).__name__}'
            ) from e


class _SimpleCData(_CData, metaclass=PyCSimpleType):

    def __init__(self, value=None):
        self._buffer_      = _alloc(self._size_)
        self._offset_      = 0
        self._b_needsfree_ = True
        if value is not None:
            self.value = value

    # ------------------------------------------------------------------
    @property
    def value(self):
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
            return _read_cstring(addr) if addr != 0 else None
        if code == 'Z':
            addr = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            return _read_cwstring(addr) if addr != 0 else None
        if code == 'O':
            idx = self._buffer_.unpack_from(_PTR_FMT, off)[0]
            return _py_object_store.get(idx)

        return self._buffer_.unpack_from(fmt, off)[0]

    @value.setter
    def value(self, v):
        code = self._type_
        fmt  = self._fmt_
        off  = self._offset_

        if code == 'c':
            if isinstance(v, int):
                v = bytes([v & 0xFF])
            if not isinstance(v, (bytes, bytearray)) or len(v) != 1:
                raise TypeError('one character bytes or integer expected')
            self._buffer_.write(v, off)
            return

        if code == 'u':
            if not isinstance(v, str) or len(v) != 1:
                raise TypeError('a unicode character is required')
            self._buffer_.write(v.encode('utf-16-le'), off)
            return

        if code in ('z', 'Z', 'P'):
            if v is None:
                self._buffer_.pack_into(_PTR_FMT, off, 0)
            elif isinstance(v, int):
                self._buffer_.pack_into(_PTR_FMT, off, v)
            elif isinstance(v, (bytes, bytearray)):
                s   = bytes(v) + b'\x00'
                tmp = _alloc(len(s), s)
                self._buffer_.pack_into(_PTR_FMT, off, tmp.address)
                if self._objects is None:
                    self._objects = {}
                self._objects[tmp.address] = tmp
            elif isinstance(v, str):
                s   = v.encode() + b'\x00'
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
        except Exception as e:
            raise OverflowError(str(e)) from e

    # ------------------------------------------------------------------
    def __repr__(self):
        try:
            return f'{type(self).__name__}({self.value!r})'
        except Exception:
            return f'{type(self).__name__}(<e>)'

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

# ---------------------------------------------------------------------------
# Field layout helpers
# ---------------------------------------------------------------------------

class _FieldDescriptor:
    """Data descriptor for a single Structure/Union field."""
    __slots__ = ('name', 'ctype', 'offset', 'size', 'bit_offset', 'bit_size')

    def __init__(self, name, ctype, offset, size, bit_offset=0, bit_size=0):
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
    """Compute field offsets; return (descriptors_dict, total_size, max_align)."""
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
            descriptors[fname] = _FieldDescriptor(fname, ctype, field_off, sz)
            if not is_union:
                offset = field_off + sz

    if is_union:
        total = max((sizeof(s[1]) for s in fields), default=0)
    else:
        total = offset

    if max_align > 1 and total % max_align:
        total += max_align - (total % max_align)

    return descriptors, total, max_align

# ---------------------------------------------------------------------------
# Array
# ---------------------------------------------------------------------------

def _make_array_type(item_type, count: int):
    key = (item_type, count)
    cached = _pointer_type_cache.get(key)
    if cached is not None:
        return cached

    sz   = sizeof(item_type) * count
    aln  = alignment(item_type)
    name = f'{item_type.__name__}_Array_{count}'
    ns   = {
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
        ba   = bytes(source)
        buf  = _alloc(len(ba), ba)
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
        else:
            tmp = self._type_(value)
            self._buffer_.write(tmp._buffer_.read(0, item_sz), off)

    def _read_item(self, index: int):
        off  = self._item_offset(index)
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

    @property
    def value(self):
        if issubclass(self._type_, _SimpleCData):
            code = self._type_._type_
            if code == 'c':
                raw = self._buffer_.read(self._offset_, self._size_)
                nul = raw.find(b'\x00')
                return raw[:nul] if nul >= 0 else raw
            if code == 'u':
                raw = self._buffer_.read(self._offset_, self._size_)
                return raw.decode('utf-16-le', errors='replace').rstrip('\x00')
        raise AttributeError(f'{type(self).__name__} has no .value')

    @value.setter
    def value(self, v):
        if issubclass(self._type_, _SimpleCData):
            code = self._type_._type_
            if code == 'c':
                if isinstance(v, str):
                    v = v.encode()
                raw = (v[:self._size_ - 1] + b'\x00').ljust(self._size_, b'\x00')
                self._buffer_.write(raw, self._offset_)
                return
            if code == 'u':
                raw = v.encode('utf-16-le')
                cap = self._size_ - 2
                raw = (raw[:cap] + b'\x00\x00').ljust(self._size_, b'\x00')
                self._buffer_.write(raw, self._offset_)
                return
        raise AttributeError(f'cannot set .value on {type(self).__name__}')

    @property
    def raw(self) -> bytes:
        return self._buffer_.read(self._offset_, self._size_)

    @raw.setter
    def raw(self, data: bytes):
        if len(data) > self._size_:
            raise ValueError('raw buffer is too large')
        self._buffer_.write(data.ljust(self._size_, b'\x00'), self._offset_)

    def __repr__(self):
        return f'<{type(self).__name__} object at 0x{addressof(self):08x}>'

# ---------------------------------------------------------------------------
# Structure
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
        ba   = bytes(source)
        buf  = _alloc(len(ba), ba)
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
        return f'<{type(self).__name__} at 0x{addressof(self):08x}>'

# ---------------------------------------------------------------------------
# Union
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
        return f'<{type(self).__name__} at 0x{addressof(self):08x}>'

# ---------------------------------------------------------------------------
# POINTER / _Pointer
# ---------------------------------------------------------------------------

def POINTER(ctype):
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

    def _get_ptr_val(self) -> int:
        return self._buffer_.unpack_from(self._fmt_, self._offset_)[0]

    def _set_ptr_val(self, addr: int):
        self._buffer_.pack_into(self._fmt_, self._offset_, addr)

    @property
    def contents(self):
        addr = self._get_ptr_val()
        if addr == 0:
            raise ValueError('NULL pointer access')
        return self._type_.from_address(addr)

    @contents.setter
    def contents(self, obj):
        self._set_ptr_val(addressof(obj))
        self._objects[id(obj)] = obj

    def __getitem__(self, index: int):
        addr = self._get_ptr_val()
        if addr == 0:
            raise ValueError('NULL pointer dereference')
        item_sz = sizeof(self._type_)
        inst    = self._type_.from_address(addr + index * item_sz)
        if issubclass(self._type_, _SimpleCData):
            return inst.value
        return inst

    def __setitem__(self, index: int, value):
        addr    = self._get_ptr_val()
        item_sz = sizeof(self._type_)
        target  = addr + index * item_sz
        buf     = _buf_from_addr(addr)
        if buf is None:
            raise ValueError('pointer points outside managed heap')
        if issubclass(self._type_, _SimpleCData):
            tmp = self._type_(value)
            buf.write(tmp._buffer_.read(0, item_sz), target - buf.address)
        else:
            raise TypeError('pointer item assignment not supported for this type')

    @property
    def value(self) -> int:
        return self._get_ptr_val()

    def __bool__(self) -> bool:
        return self._get_ptr_val() != 0

    def __repr__(self):
        return (f'<{type(self).__name__} at 0x{addressof(self):08x}'
                f' -> 0x{self._get_ptr_val():08x}>')

# ---------------------------------------------------------------------------
# Virtual library / FFI layer
# ---------------------------------------------------------------------------

_PROC_HANDLE   = 0
_vlibs:          dict = {_PROC_HANDLE: {}}
_vlib_names:     dict = {_PROC_HANDLE: '<process>'}
_vcallables:     dict = {}
_vlib_counter        = 0
_vaddr_counter       = 0x7FFF_0000      # 32-bit synthetic address space

def _next_vaddr() -> int:
    global _vaddr_counter
    _vaddr_counter += 8
    return _vaddr_counter

def register_library(name: str, symbols: dict) -> int:
    """
    Register Python callables as a virtual native library.

    Symbols are also mirrored onto *window* so that JS code (or other Brython
    modules using window.xxx) can call them directly from the browser.
    """
    global _vlib_counter
    _vlib_counter += 1
    h    = _vlib_counter
    syms = {}
    for sym_name, callable_ in symbols.items():
        addr = _next_vaddr()
        syms[sym_name] = (callable_, addr)
        _vcallables[addr] = callable_
        # Brython-specific: expose on window for JS interop
        try:
            setattr(window, f'__vlib_{sym_name}', callable_)
        except Exception:
            pass
    _vlibs[h]      = syms
    _vlib_names[h] = name
    return h

def _vlib_sym(lib_handle: int, name: str) -> int:
    syms = _vlibs.get(lib_handle)
    if syms is None:
        raise OSError(f'invalid library handle {lib_handle}')
    entry = syms.get(name)
    if entry is None:
        if lib_handle == _PROC_HANDLE:
            addr = _next_vaddr()
            syms[name] = (None, addr)
            return addr
        raise OSError(f'symbol {name!r} not found in '
                      f'{_vlib_names[lib_handle]!r}')
    return entry[1]

def dlopen(name=None, mode: int = RTLD_LOCAL) -> int:
    if name is None:
        return _PROC_HANDLE
    for h, n in _vlib_names.items():
        if n == name or _os.path.basename(n) == _os.path.basename(name):
            return h
    raise OSError(f'dlopen({name!r}): no registered implementation — '
                  f'use register_library() to provide one')

def dlclose(handle: int) -> None:
    _vlibs.pop(handle, None)
    _vlib_names.pop(handle, None)

def dlsym(handle: int, name: str) -> int:
    return _vlib_sym(handle, name)

# ---------------------------------------------------------------------------
# CFuncPtr
# ---------------------------------------------------------------------------

class PyCFuncPtrType(type):
    def __mul__(cls, count: int):
        return _make_array_type(cls, count)


class CFuncPtr(_CData, metaclass=PyCFuncPtrType):
    """Base class for ctypes function pointer types."""

    _size_   = _PTR_SZ
    _align_  = _PTR_SZ
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
            self._callable_ = first
            addr = _next_vaddr()
            _vcallables[addr] = first
            self._func_addr_  = addr
            self._buffer_.pack_into(_PTR_FMT, 0, addr)
            # Brython: expose on window so JS can trigger callbacks
            try:
                setattr(window, f'__cb_{addr}', first)
            except Exception:
                pass

        elif isinstance(first, int):
            self._func_addr_ = first
            self._buffer_.pack_into(_PTR_FMT, 0, first)

        elif isinstance(first, tuple) and len(first) == 2:
            name_or_ord, lib = first
            if isinstance(name_or_ord, str):
                lib_handle = lib if isinstance(lib, int) \
                             else getattr(lib, '_handle', id(lib))
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
                f'no Python callable at 0x{addr:08x} — '
                f'use register_library() / CFuncPtr(callback)')

        py_args = []
        for arg in args:
            if isinstance(arg, _CArgObject):
                py_args.append(arg._obj)
            elif isinstance(arg, _SimpleCData):
                py_args.append(arg.value)
            else:
                py_args.append(arg)

        result = cb(*py_args)

        if self.errcheck is not None:
            result = self.errcheck(result, self, args)

        if self.restype is not None and not isinstance(result, self.restype):
            try:
                result = self.restype(result).value
            except Exception:
                pass

        return result

    def __repr__(self):
        return f'<{type(self).__name__} at 0x{addressof(self):08x}>'

# ---------------------------------------------------------------------------
# Low-level call functions
# ---------------------------------------------------------------------------

def call_function(func_addr: int, arguments: tuple):
    cb = _vcallables.get(func_addr)
    if cb is None:
        raise OSError(f'call_function: no callable at 0x{func_addr:08x}')
    return cb(*arguments)

def call_cdeclfunction(func_addr: int, arguments: tuple):
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
    return _py_object_store.get(addr)

def Py_INCREF(obj): pass
def Py_DECREF(obj): pass

# ---------------------------------------------------------------------------
# _impl_* helpers  +  their synthetic addresses
# ---------------------------------------------------------------------------

def _impl_cast(obj, obj2, typ):
    return cast(obj2 if obj2 is not None else obj, typ)

def _impl_memmove(dst, src, count):
    def _addr(x):
        if isinstance(x, int):          return x
        if isinstance(x, _CData):       return addressof(x)
        if isinstance(x, _CArgObject):  return addressof(x._obj) + x._offset
        raise TypeError(f'memmove: bad argument {type(x)}')
    dst_addr = _addr(dst)
    src_addr = _addr(src)
    n        = int(count)
    src_buf  = _buf_from_addr(src_addr)
    dst_buf  = _buf_from_addr(dst_addr)
    data     = src_buf.read(src_addr - src_buf.address, n) if src_buf else bytes(n)
    if dst_buf:
        dst_buf.write(data, dst_addr - dst_buf.address)
    return dst_addr

def _impl_memset(dst, c, count):
    def _addr(x):
        if isinstance(x, int):    return x
        if isinstance(x, _CData): return addressof(x)
        raise TypeError(f'memset: bad argument {type(x)}')
    dst_addr = _addr(dst)
    n        = int(count)
    val      = int(c) & 0xFF
    dst_buf  = _buf_from_addr(dst_addr)
    if dst_buf:
        dst_buf.write(bytes([val] * n), dst_addr - dst_buf.address)
    return dst_addr

def _impl_string_at(addr, size=-1):
    if isinstance(addr, _CData):
        addr = addressof(addr)
    buf = _buf_from_addr(int(addr))
    if buf is None:
        return b''
    off = int(addr) - buf.address
    raw = buf.read(off)
    if size < 0:
        nul = raw.find(b'\x00')
        return raw[:nul] if nul >= 0 else raw
    return raw[:size]

def _impl_wstring_at(addr, size=-1):
    if isinstance(addr, _CData):
        addr = addressof(addr)
    wsz = 2
    enc = 'utf-16-le'
    buf = _buf_from_addr(int(addr))
    if buf is None:
        return ''
    off = int(addr) - buf.address
    raw = buf.read(off)
    if size < 0:
        i = 0
        while i + wsz <= len(raw):
            if raw[i:i+wsz] == b'\x00\x00':
                break
            i += wsz
        raw = raw[:i]
    else:
        raw = raw[:size * wsz]
    return raw.decode(enc, errors='replace')

# Allocate stable synthetic addresses and wire up callables
_cast_addr       = _next_vaddr()
_memmove_addr    = _next_vaddr()
_memset_addr     = _next_vaddr()
_string_at_addr  = _next_vaddr()
_wstring_at_addr = _next_vaddr()

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
# Brython-specific JS bridge helpers
# ---------------------------------------------------------------------------

def js_array_buffer(obj) -> object:
    """
    Wrap a ctypes buffer (or any _CData instance) as a JavaScript
    ArrayBuffer so that Web APIs (WebGL, WebAudio, fetch body, etc.)
    can consume it without copying.

    Returns a JS ArrayBuffer whose underlying bytes mirror obj's _RawBuf.
    Changes to the returned ArrayBuffer are NOT reflected back into the
    Python-side buffer — use js_data_view for bidirectional access.
    """
    if isinstance(obj, _CData):
        raw = obj._buffer_.read(obj._offset_, sizeof(obj))
    elif isinstance(obj, (bytes, bytearray)):
        raw = bytes(obj)
    else:
        raise TypeError(f'js_array_buffer: cannot wrap {type(obj)}')

    # Build a Uint8Array from individual bytes, then return its .buffer
    uint8 = window.Uint8Array.new(len(raw))
    for i, b in enumerate(raw):
        uint8[i] = b
    return uint8.buffer


def js_data_view(obj) -> object:
    """
    Return a JS DataView backed by a fresh ArrayBuffer that is
    pre-loaded with obj's current bytes.

    Write-back example::

        dv = js_data_view(my_struct)
        dv.setInt32(0, 42, True)     # True = little-endian
        # then sync back if needed:  js_sync_back(dv, my_struct)
    """
    ab = js_array_buffer(obj)
    return window.DataView.new(ab)


def js_sync_back(dv, obj):
    """
    Copy bytes from a JS DataView (or Uint8Array / ArrayBuffer) back
    into a ctypes object's _RawBuf.

    Use after modifying a DataView returned by js_data_view().
    """
    n = sizeof(obj) if isinstance(obj, _CData) else len(obj._buffer_)
    # dv.buffer is the underlying ArrayBuffer; wrap it as Uint8Array to iterate
    try:
        src = window.Uint8Array.new(dv.buffer)
    except Exception:
        src = window.Uint8Array.new(dv)   # maybe dv already is a Uint8Array

    buf  = obj._buffer_ if isinstance(obj, _CData) else obj
    off  = obj._offset_ if isinstance(obj, _CData) else 0
    data = bytes(src[i] for i in range(n))
    buf.write(data, off)


def js_callback(func, restype=None, argtypes=None):
    """
    Wrap a Python callable as a JavaScript-callable function and
    register it in the virtual callable table.

    Returns a CFuncPtr whose .value is the synthetic address, and also
    exposes the callback on window.__pycb_<addr> for direct JS use.

    Example::

        @js_callback
        def on_message(ptr, length):
            data = _impl_string_at(ptr, length)
            print('got:', data)

        window.myWebSocket.onmessage = window.__pycb_<addr>
    """
    fptr = CFuncPtr(func)
    if restype is not None:
        fptr.restype = restype
    if argtypes is not None:
        fptr.argtypes = argtypes
    addr = fptr._func_addr_
    try:
        setattr(window, f'__pycb_{addr}', func)
    except Exception:
        pass
    return fptr

# ---------------------------------------------------------------------------
# Public API surface
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
    # Brython-specific JS bridge
    'js_array_buffer', 'js_data_view', 'js_sync_back', 'js_callback',
]

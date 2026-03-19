from __future__ import annotations
"""
_numpy_binary_shims.py
Pure-Python shims for every compiled (.so/.pyd) extension module that
NumPy 2.x ships.  Drop this file next to (or import before) numpy so that
the normal `import numpy` succeeds in environments like Brython where
no C extensions can be loaded.
"""
import sys
import math
import cmath
import operator
import random as _random
import struct  as _struct
import array   as _array
import itertools
import functools
import weakref
import types
import re
import abc
from collections import namedtuple
import collections.abc as _cabc

# Python 3.8 compatibility: NumPy's pure-Python code may reference
# `types.GenericAlias`. Provide a minimal stub when missing.
if not hasattr(types, "GenericAlias"):
    class GenericAlias:  # pragma: no cover
        def __init__(self, origin, params):
            self.__origin__ = origin
            self.__params__ = params

        def __repr__(self):
            return f"{self.__origin__}[{self.__params__!r}]"

    types.GenericAlias = GenericAlias

# ---------------------------------------------------------------------------
# Helper: build a fake module and register it in sys.modules
# ---------------------------------------------------------------------------
def _make_module(fullname: str, **attrs) -> types.ModuleType:
    mod = sys.modules.get(fullname)
    if mod is None:
        mod = types.ModuleType(fullname)
        mod.__name__     = fullname
        mod.__package__  = '.'.join(fullname.split('.')[:-1])
        mod.__loader__   = None
        mod.__spec__     = None
        mod.__file__     = '<_numpy_binary_shims>'
        for k, v in attrs.items():
            setattr(mod, k, v)
        sys.modules[fullname] = mod
    parts = fullname.split('.')
    if len(parts) > 1:
        parent = sys.modules.get('.'.join(parts[:-1]))
        if parent is not None:
            try:
                parent.__dict__[parts[-1]] = mod
            except (AttributeError, TypeError):
                pass
    return mod

_aliases_injected = False
_ensure_running   = False

def _ensure_submodule_attrs():
    global _aliases_injected, _ensure_running
    if _ensure_running:
        return
    _ensure_running = True
    try:
        _ensure_submodule_attrs_inner()
    finally:
        _ensure_running = False

def _ensure_submodule_attrs_inner():
    global _aliases_injected
    for fullname, mod in list(sys.modules.items()):
        if not fullname.startswith('numpy'):
            continue
        parts = fullname.split('.')
        if len(parts) < 2:
            continue
        parent = sys.modules.get('.'.join(parts[:-1]))
        if parent is not None:
            try:
                if parts[-1] not in parent.__dict__:
                    parent.__dict__[parts[-1]] = mod
            except (AttributeError, TypeError):
                pass
    if not _aliases_injected:
        _np = sys.modules.get('numpy')
        if _np is not None:
            _nd = _np.__dict__
            _builtin_aliases = {
                'bool': bool, 'int': int, 'float': float,
                'complex': complex, 'object': object, 'str': str,
                'bytes': bytes, 'unicode': str,
            }
            _scalar_aliases = {
                'cfloat': complex64, 'cdouble': complex128, 'csingle': complex64,
                'clongdouble': complex128, 'longcomplex': complex128,
                'float_': float64, 'int_': int64, 'complex_': complex128,
                'bool8': bool_, 'int0': int64, 'uint0': uint64,
                'string_': bytes_, 'unicode_': str_, 'object_': object_,
                'longdouble': float64, 'long': int64, 'ulong': uint64,
            }
            for alias, val in {**_builtin_aliases, **_scalar_aliases}.items():
                if alias not in _nd:
                    _nd[alias] = val
            _aliases_injected = True

# ===========================================================================
# Section 1 – dtype
# ===========================================================================
_DTYPE_TABLE = {
    'b': ('b', 1, int), 'B': ('B', 1, int),
    'h': ('h', 2, int), 'H': ('H', 2, int),
    'i': ('i', 4, int), 'I': ('I', 4, int),
    'l': ('l', 8, int), 'L': ('L', 8, int),
    'q': ('q', 8, int), 'Q': ('Q', 8, int),
    'e': ('e', 2, float), 'f': ('f', 4, float),
    'd': ('d', 8, float), 'g': ('d', 8, float),
    'F': ('f', 4, complex), 'D': ('d', 8, complex),
    'G': ('d', 8, complex), '?': ('?', 1, bool),
    'S': ('s', 1, bytes), 'U': ('U', 4, str),
    'O': ('P', 8, object), 'V': ('B', 1, bytes),
    'M': ('q', 8, int), 'm': ('q', 8, int),
}

STR_TO_KIND = {
    'bool': '?', 'int8': 'b', 'byte': 'b',
    'uint8': 'B', 'ubyte': 'B',
    'int16': 'h', 'short': 'h',
    'uint16': 'H', 'ushort': 'H',
    'int32': 'i', 'intc': 'i',
    'uint32': 'I', 'uintc': 'I',
    'int64': 'q', 'int_': 'q', 'long': 'q', 'intp': 'q',
    'uint64': 'Q', 'uint': 'Q', 'ulong': 'Q', 'uintp': 'Q',
    'longlong': 'q', 'ulonglong': 'Q',
    'float16': 'e', 'half': 'e',
    'float32': 'f', 'single': 'f',
    'float64': 'd', 'double': 'd', 'float_': 'd', 'float': 'd',
    'longdouble': 'g',
    'complex64': 'F',
    'complex128': 'D', 'complex_': 'D', 'cdouble': 'D', 'complex': 'D',
    'clongdouble': 'G',
    'bytes_': 'S', 'string_': 'S',
    'str_': 'U', 'unicode_': 'U',
    'object_': 'O', 'object': 'O',
    'void': 'V',
    'datetime64': 'M', 'timedelta64': 'm',
    'b1': '?', 'i1': 'b', 'u1': 'B', 'i2': 'h', 'u2': 'H',
    'i4': 'i', 'u4': 'I', 'i8': 'q', 'u8': 'Q',
    'f2': 'e', 'f4': 'f', 'f8': 'd',
    'c8': 'F', 'c16': 'D',
}

_ITEMSIZES = {
    '?': 1, 'b': 1, 'B': 1, 'h': 2, 'H': 2, 'e': 2,
    'i': 4, 'I': 4, 'f': 4, 'F': 8,
    'q': 8, 'Q': 8, 'd': 8, 'D': 16,
    'g': 8, 'G': 16, 'l': 8, 'L': 8,
    'S': 1, 'U': 4, 'O': 8, 'V': 1, 'M': 8, 'm': 8,
}

KIND_NAMES = {
    '?': 'bool', 'b': 'int8', 'B': 'uint8', 'h': 'int16', 'H': 'uint16',
    'i': 'int32', 'I': 'uint32', 'q': 'int64', 'Q': 'uint64',
    'e': 'float16', 'f': 'float32', 'd': 'float64', 'g': 'float64',
    'F': 'complex64', 'D': 'complex128', 'G': 'complex128',
    'S': 'bytes', 'U': 'str_', 'O': 'object_', 'V': 'void',
    'M': 'datetime64', 'm': 'timedelta64',
}

_ABSTRACT_KIND = {
    '?': 'b', 'b': 'i', 'h': 'i', 'i': 'i', 'l': 'i', 'q': 'i',
    'B': 'u', 'H': 'u', 'I': 'u', 'L': 'u', 'Q': 'u',
    'e': 'f', 'f': 'f', 'd': 'f', 'g': 'f',
    'F': 'c', 'D': 'c', 'G': 'c',
    'S': 'S', 'U': 'U', 'O': 'O', 'V': 'V', 'M': 'M', 'm': 'm',
}

_NUMPY_CHAR = {
    '?': '?', 'b': 'b', 'h': 'h', 'i': 'i', 'l': 'l', 'q': 'q',
    'B': 'B', 'H': 'H', 'I': 'I', 'L': 'L', 'Q': 'Q',
    'e': 'e', 'f': 'f', 'd': 'd', 'g': 'g',
    'F': 'F', 'D': 'D', 'G': 'G',
    'S': 'S', 'U': 'U', 'O': 'O', 'V': 'V', 'M': 'M', 'm': 'm',
}

class dtype:
    _cache: dict = {}
    alignment = 1
    byteorder = '='
    descr = []
    fields = None
    flags = 0
    hasobject = False
    isalignedstruct = False
    isbuiltin = 1
    isnative = True
    itemsize = 8
    kind = 'f'
    metadata = None
    name = 'float64'
    names = None
    ndim = 0
    num = 12
    shape = ()
    str = '<f8'
    subdtype = None
    type = float

    def __new__(cls, spec, align=False, copy=False):
        if isinstance(spec, dtype):
            return spec
        if isinstance(spec, type):
            _map = {bool: '?', int: 'q', float: 'd', complex: 'D',
                    bytes: 'S', str: 'U', object: 'O'}
            if spec in _map:
                spec = _map[spec]
            else:
                _sc = [
                    (bool_, '?'), (datetime64, 'M'), (timedelta64, 'm'),
                    (float16, 'e'), (float32, 'f'), (float64, 'd'),
                    (complex64, 'F'), (complex128, 'D'),
                    (int8, 'b'), (uint8, 'B'), (int16, 'h'), (uint16, 'H'),
                    (int32, 'i'), (uint32, 'I'), (int64, 'q'), (uint64, 'Q'),
                    (bytes_, 'S'), (str_, 'U'), (void, 'V'), (object_, 'O'),
                    (integer, 'q'), (floating, 'd'), (complexfloating, 'D'),
                    (generic, 'O'),
                ]
                # Determine the dtype "kind" character from the passed scalar type.
                # Important: do not overwrite `spec` before the `issubclass()` checks,
                # otherwise the mapping loop can never match.
                _kind_char = 'O'
                for _sc_cls, _sc_char in _sc:
                    try:
                        if issubclass(spec, _sc_cls):
                            _kind_char = _sc_char
                            break
                    except TypeError:
                        pass
                spec = _kind_char
        if not isinstance(spec, str):
            spec = str(spec)
        byteorder = '='
        if spec and spec[0] in '<>=|':
            byteorder = spec[0]
            spec = spec[1:]
        extra = 0
        m = re.match(r'^([A-Za-z_]\w*?)(\d+)$', spec)
        if m:
            base, num = m.group(1), int(m.group(2))
            if base in ('S', 'U', 'V', 'a'):
                extra = num
                spec = base
            elif base in STR_TO_KIND:
                spec = STR_TO_KIND.get(base + str(num), STR_TO_KIND.get(base, 'O'))
        elif spec in STR_TO_KIND:
            spec = STR_TO_KIND[spec]
        kind = spec if len(spec) == 1 else STR_TO_KIND.get(spec, 'O')
        key = (kind, byteorder, extra)
        cached = cls._cache.get(key)
        if cached is not None and not copy:
            return cached
        obj = object.__new__(cls)
        obj._char = _NUMPY_CHAR.get(kind, kind)
        obj.kind = _ABSTRACT_KIND.get(kind, kind)
        obj.byteorder = byteorder
        obj._extra = extra
        base_sz = _ITEMSIZES.get(kind, 1)
        obj.itemsize = (extra or base_sz) * (4 if kind == 'U' else 1)
        obj.name = (KIND_NAMES.get(kind, 'object') + (str(extra) if extra else ''))
        obj.str = byteorder + kind + str(obj.itemsize)
        obj.num = list(_ITEMSIZES.keys()).index(kind) if kind in _ITEMSIZES else 0
        obj.fields = None
        obj.names = None
        obj.subdtype = None
        obj.shape = ()
        obj.ndim = 0
        obj.isbuiltin = 1
        obj.isnative = True
        obj.descr = [(obj.name, obj.str)]
        obj.alignment = obj.itemsize
        obj.flags = 0
        obj.hasobject = (kind == 'O')
        obj.type = None
        cls._cache[key] = obj
        return obj

    @property
    def char(self):
        return getattr(self, "_char", self.kind)

    def __eq__(self, other):
        if isinstance(other, dtype):
            return self.str == other.str
        try:
            return self == dtype(other)
        except Exception:
            return False

    def __hash__(self):
        return hash(self.str)

    def __repr__(self):
        return f"dtype('{self.name}')"

    def __str__(self):
        return self.name

    def newbyteorder(self, new_order='S'):
        return dtype(new_order + self.kind)

    @staticmethod
    def _from_typestr(ts: str) -> 'dtype':
        return dtype(ts)

    @classmethod
    def from_fields(cls, fields_dict: dict) -> 'dtype':
        obj = object.__new__(cls)
        obj.kind = 'V'
        obj.byteorder = '|'
        obj._extra = 0
        offset = 0
        field_map = {}
        names = []
        for fname, fdt in fields_dict.items():
            fdt = dtype(fdt)
            field_map[fname] = (fdt, offset)
            names.append(fname)
            offset += fdt.itemsize
        obj.itemsize = offset
        obj.fields = field_map
        obj.names = tuple(names)
        obj.name = 'void' + str(offset * 8)
        obj.str = '|V' + str(offset)
        obj.num = 0
        obj.subdtype = None
        obj.shape = ()
        obj.ndim = 0
        obj.isbuiltin = 0
        obj.isnative = True
        obj.descr = [(n, str(dtype(v))) for n, v in fields_dict.items()]
        obj.alignment = 1
        obj.flags = 0
        obj.hasobject = False
        obj.type = object
        return obj

    @property
    def base(self):
        if self.subdtype:
            return self.subdtype[0]
        return self

    @classmethod
    def __class_getitem__(cls, item):
        return types.GenericAlias(cls, (item,) if not isinstance(item, tuple) else item)

    def __lt__(self, other): return NotImplemented
    def __le__(self, other): return NotImplemented
    def __gt__(self, other): return NotImplemented
    def __ge__(self, other): return NotImplemented

_CHAR_TO_SCALAR = None

def _dtype_patch_types():
    global _CHAR_TO_SCALAR
    _CHAR_TO_SCALAR = {
        '?': bool_,
        'b': int8, 'B': uint8,
        'h': int16, 'H': uint16,
        'i': int32, 'I': uint32,
        'q': int64, 'Q': uint64,
        'l': int64, 'L': uint64,
        'e': float16, 'f': float32, 'd': float64, 'g': float64,
        'F': complex64, 'D': complex128, 'G': complex128,
        'S': bytes_, 'U': str_, 'O': object_, 'V': void,
        'M': datetime64, 'm': timedelta64,
    }
    for key, obj in dtype._cache.items():
        kind = key[0]
        obj.type = _CHAR_TO_SCALAR.get(kind, object_)

bool_dtype = dtype('?')
int8_dtype = dtype('b')
int16_dtype = dtype('h')
int32_dtype = dtype('i')
int64_dtype = dtype('q')
uint8_dtype = dtype('B')
uint16_dtype = dtype('H')
uint32_dtype = dtype('I')
uint64_dtype = dtype('Q')
float16_dtype = dtype('e')
float32_dtype = dtype('f')
float64_dtype = dtype('d')
complex64_dtype = dtype('F')
complex128_dtype = dtype('D')
object_dtype = dtype('O')

# ===========================================================================
# Section 2 – Scalar type hierarchy
# ===========================================================================
class generic:
    type = None
    kind = ''
    def __init__(self, value=0):
        self._value = self._coerce(value)
    def _coerce(self, v):
        return v
    @property
    def flat(self): return iter([self._value])
    def item(self): return self._value
    def tolist(self): return self._value
    def __repr__(self): return f'{type(self).__name__}({self._value!r})'
    def __str__(self): return str(self._value)
    def __eq__(self, o): return self._value == (o._value if isinstance(o, generic) else o)
    def __hash__(self): return hash(self._value)
    def __add__(self, o): return self._value + (o._value if isinstance(o, generic) else o)
    def __radd__(self, o): return (o._value if isinstance(o, generic) else o) + self._value
    def __mul__(self, o): return self._value * (o._value if isinstance(o, generic) else o)
    def __rmul__(self, o): return (o._value if isinstance(o, generic) else o) * self._value
    def __sub__(self, o): return self._value - (o._value if isinstance(o, generic) else o)
    def __rsub__(self, o): return (o._value if isinstance(o, generic) else o) - self._value
    def __truediv__(self, o): return self._value / (o._value if isinstance(o, generic) else o)
    def __floordiv__(self, o): return self._value // (o._value if isinstance(o, generic) else o)
    def __mod__(self, o): return self._value % (o._value if isinstance(o, generic) else o)
    def __pow__(self, o): return self._value ** (o._value if isinstance(o, generic) else o)
    def __neg__(self): return -self._value
    def __pos__(self): return +self._value
    def __abs__(self): return abs(self._value)
    def __bool__(self): return bool(self._value)
    def __int__(self): return int(self._value)
    def __float__(self): return float(self._value)
    def __complex__(self): return complex(self._value)
    def __lt__(self, o): return self._value < (o._value if isinstance(o, generic) else o)
    def __le__(self, o): return self._value <= (o._value if isinstance(o, generic) else o)
    def __gt__(self, o): return self._value > (o._value if isinstance(o, generic) else o)
    def __ge__(self, o): return self._value >= (o._value if isinstance(o, generic) else o)
    @property
    def dtype(self): return dtype('d')
    @property
    def shape(self): return ()
    @property
    def ndim(self): return 0
    @property
    def size(self): return 1
    @property
    def T(self): return self
    def astype(self, dt): return type(self)(self._value)
    def reshape(self, *args): return self
    def copy(self, order='C'): return type(self)(self._value)
    @property
    def data(self): return memoryview(bytes([int(self._value) & 0xFF]))
    @property
    def flags(self): return {}
    @property
    def itemsize(self): return 8
    @property
    def strides(self): return ()
    @property
    def base(self): return None
    @property
    def real(self): return self._value.real if hasattr(self._value, 'real') else self._value
    @property
    def imag(self): return self._value.imag if hasattr(self._value, 'imag') else 0
    @property
    def nbytes(self): return 8
    def __array_namespace__(self, *, api_version=None): return None
    def __copy__(self): return type(self)(self._value)
    def __deepcopy__(self, memo=None): return type(self)(self._value)
    def all(self, axis=None, out=None, keepdims=False, **kw): return bool(self._value)
    def any(self, axis=None, out=None, keepdims=False, **kw): return bool(self._value)
    def argmax(self, axis=None, out=None, **kw): return 0
    def argmin(self, axis=None, out=None, **kw): return 0
    def argsort(self, axis=-1, kind=None, order=None, **kw): return 0
    def byteswap(self, inplace=False): return type(self)(self._value)
    def choose(self, choices, out=None, mode='raise'): return choices[int(self._value)]
    def clip(self, min=None, max=None, out=None, **kw):
        v = self._value
        if min is not None and v < min: v = min
        if max is not None and v > max: v = max
        return type(self)(v)
    def compress(self, condition, axis=None, out=None): return self
    def conj(self): return type(self)(self._value)
    def conjugate(self): return type(self)(self._value)
    def cumprod(self, axis=None, dtype=None, out=None): return type(self)(self._value)
    def cumsum(self, axis=None, dtype=None, out=None): return type(self)(self._value)
    def diagonal(self, offset=0, axis1=0, axis2=1): return self
    def dump(self, file):
        import pickle
        with open(file, 'wb') as f: pickle.dump(self, f)
    def dumps(self):
        import pickle
        return pickle.dumps(self)
    def fill(self, value): self._value = self._coerce(value)
    def flatten(self, order='C'): return self
    def getfield(self, dtype, offset=0): return self
    def setfield(self, val, dtype, offset=0): pass
    def setflags(self, write=None, align=None, uic=None): pass
    def max(self, axis=None, out=None, **kw): return self._value
    def mean(self, axis=None, dtype=None, out=None, **kw): return float(self._value)
    def min(self, axis=None, out=None, **kw): return self._value
    def nonzero(self): return ([0],) if self._value else ([],)
    def prod(self, axis=None, dtype=None, out=None, **kw): return self._value
    def put(self, indices, values, mode='raise'): pass
    def ravel(self, order='C'): return self
    def repeat(self, repeats, axis=None): return self
    def resize(self, *new_shape, refcheck=True): pass
    def round(self, decimals=0, out=None): return type(self)(round(float(self._value), decimals))
    def searchsorted(self, v, side='left', sorter=None): return 0
    def sort(self, axis=-1, kind=None, order=None, **kw): pass
    def squeeze(self, axis=None): return self
    def std(self, axis=None, dtype=None, out=None, ddof=0, **kw): return 0.0
    def as_integer_ratio(self): return float(self._value).as_integer_ratio()
    def is_integer(self): return float(self._value).is_integer()
    def bit_count(self): return bin(int(self._value)).count("1")
    def sum(self, axis=None, dtype=None, out=None, **kw): return self._value
    def swapaxes(self, axis1, axis2): return self
    def take(self, indices, axis=None, out=None, mode='raise'): return self
    def to_device(self, device, *, stream=None): return self
    def tobytes(self, order='C'): return bytes([int(self._value) & 0xFF])
    def tofile(self, fid, sep='', format='%s'): pass
    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None): return self._value
    def transpose(self, *axes): return self
    def var(self, axis=None, dtype=None, out=None, ddof=0, **kw): return 0.0
    def view(self, dtype=None, type=None): return self

class number(generic):
    type = None
    kind = ''
    @classmethod
    def __class_getitem__(cls, item):
        return types.GenericAlias(cls, (item,) if not isinstance(item, tuple) else item)

class integer(number):
    type = None
    kind = ''
    def is_integer(self): return True
    @classmethod
    def __class_getitem__(cls, item):
        return types.GenericAlias(cls, (item,) if not isinstance(item, tuple) else item)

class signedinteger(integer):
    type = None
    kind = ''

class unsignedinteger(integer):
    type = None
    kind = ''

class inexact(number):
    type = None
    kind = ''

class floating(inexact):
    type = None
    kind = ''

class complexfloating(inexact):
    type = None
    kind = ''

class flexible(generic):
    type = None
    kind = ''

class character(flexible):
    type = None
    kind = ''

generic.type = generic
number.type = number
integer.type = integer
signedinteger.type = signedinteger
unsignedinteger.type = unsignedinteger
inexact.type = inexact
floating.type = floating
complexfloating.type = complexfloating
flexible.type = flexible
character.type = character

def _make_int_type(name, kind, bits, signed=True):
    lo = -(1 << (bits - 1)) if signed else 0
    hi = (1 << (bits - 1)) - 1 if signed else (1 << bits) - 1
    base = signedinteger if signed else unsignedinteger
    def _coerce(self, v):
        v = int(v)
        return lo + (v - lo) % (1 << bits)
    cls = type(name, (base,), {
        '_coerce': _coerce, '_kind': kind, '_bits': bits,
        '_lo': lo, '_hi': hi,
        'dtype': property(lambda self: dtype(kind)),
    })
    return cls

def _make_float_type(name, kind):
    base = floating
    cls = type(name, (base,), {
        '_coerce': lambda self, v: float(v),
        '_kind': kind,
        'dtype': property(lambda self: dtype(kind)),
    })
    return cls

def _make_complex_type(name, kind):
    base = complexfloating
    cls = type(name, (base,), {
        '_coerce': lambda self, v: complex(v),
        '_kind': kind,
        'dtype': property(lambda self: dtype(kind)),
    })
    return cls

bool_ = _make_int_type('bool_', '?', 1, signed=False)
int8 = _make_int_type('int8', 'b', 8, signed=True)
int16 = _make_int_type('int16', 'h', 16, signed=True)
int32 = _make_int_type('int32', 'i', 32, signed=True)
int64 = _make_int_type('int64', 'q', 64, signed=True)
uint8 = _make_int_type('uint8', 'B', 8, signed=False)
uint16 = _make_int_type('uint16', 'H', 16, signed=False)
uint32 = _make_int_type('uint32', 'I', 32, signed=False)
uint64 = _make_int_type('uint64', 'Q', 64, signed=False)
int_ = int64
intp = int64
intc = int32
byte = int8
short = int16
long = int64
longlong = int64
ubyte = uint8
ushort = uint16
uintc = uint32
ulong = uint64
ulonglong = uint64
uintp = uint64
float16 = _make_float_type('float16', 'e')
float32 = _make_float_type('float32', 'f')
float64 = _make_float_type('float64', 'd')
float_ = float64
double = float64
longdouble = float64
single = float32
half = float16
complex64 = _make_complex_type('complex64', 'F')
complex128 = _make_complex_type('complex128', 'D')
complex_ = complex128
cdouble = complex128
clongdouble = complex128
bytes_ = type('bytes_', (character,), {'_coerce': lambda self, v: bytes(v) if not isinstance(v, bytes) else v})
str_ = type('str_', (character,), {'_coerce': lambda self, v: str(v)})
object_ = type('object_', (generic,), {'_coerce': lambda self, v: v})
void = type('void', (flexible,), {'_coerce': lambda self, v: bytes(v) if not isinstance(v, bytes) else v})

class datetime64(generic):
    def __init__(self, value=0, unit=None):
        self._value = value
        self._unit = unit
    def _coerce(self, v): return int(v) if not isinstance(v, str) else v
    def __str__(self): return str(self._value)
    def __repr__(self): return f'numpy.datetime64({self._value!r})'

class timedelta64(generic):
    def __init__(self, value=0, unit=None):
        self._value = value
        self._unit = unit
    def _coerce(self, v): return int(v)
    def __str__(self): return str(self._value)
    def __repr__(self): return f'numpy.timedelta64({self._value!r})'

_dtype_patch_types()

# ===========================================================================
# Section 3 – ndarray
# ===========================================================================
def _prod(seq):
    r = 1
    for x in seq:
        r *= x
    return r

def _flat_list(nested):
    if isinstance(nested, (list, tuple)):
        out = []
        for x in nested:
            out.extend(_flat_list(x))
        return out
    return [nested]

def _nested_shape(obj):
    if not isinstance(obj, (list, tuple)):
        return ()
    if len(obj) == 0:
        return (0,)
    inner = _nested_shape(obj[0])
    return (len(obj),) + inner

def _coerce_dtype(dt):
    """Convert various type specifications to a dtype instance."""
    if dt is None:
        return None
    if isinstance(dt, dtype):
        return dt

    # --- FIX: Direct type name matching for dynamically created scalar types ---
    # This handles cases where issubclass() fails with dynamic types created via type()
    _name_to_char = {
        'bool_': '?', 'int8': 'b', 'uint8': 'B',
        'int16': 'h', 'uint16': 'H',
        'int32': 'i', 'uint32': 'I',
        'int64': 'q', 'uint64': 'Q',
        'float16': 'e', 'float32': 'f', 'float64': 'd',
        'complex64': 'F', 'complex128': 'D',
        'bytes_': 'S', 'str_': 'U', 'object_': 'O',
        'void': 'V', 'datetime64': 'M', 'timedelta64': 'm',
        # Aliases
        'int_': 'q', 'intp': 'q', 'intc': 'i',
        'uint': 'Q', 'uintp': 'Q', 'uintc': 'I',
        'float_': 'd', 'double': 'd', 'single': 'f', 'half': 'e',
        'longdouble': 'g', 'complex_': 'D', 'cdouble': 'D',
        'byte': 'b', 'short': 'h', 'long': 'q', 'longlong': 'q',
        'ubyte': 'B', 'ushort': 'H', 'ulong': 'Q', 'ulonglong': 'Q',
    }

    # Check by type __name__ attribute (works for dynamically created types)
    if hasattr(dt, '__name__'):
        name = dt.__name__
        if name in _name_to_char:
            return dtype(_name_to_char[name])

    if isinstance(dt, type):
        # Fast path: Python builtins
        _builtin_map = {bool: '?', int: 'q', float: 'd', complex: 'D',
                        bytes: 'S', str: 'U', object: 'O'}
        if dt in _builtin_map:
            return dtype(_builtin_map[dt])

        # Fallback: try issubclass (may fail with dynamic types)
        _hierarchy = [
            (bool_,      '?'),
            (datetime64, 'M'), (timedelta64, 'm'),
            (float16,    'e'), (float32,    'f'), (float64,    'd'),
            (complex64,  'F'), (complex128, 'D'),
            (int8,  'b'), (uint8,  'B'),
            (int16, 'h'), (uint16, 'H'),
            (int32, 'i'), (uint32, 'I'),
            (int64, 'q'), (uint64, 'Q'),
            (bytes_, 'S'), (str_, 'U'), (void, 'V'), (object_, 'O'),
            (integer,     'q'), (floating,    'd'), (complexfloating, 'D'),
            (generic,     'O'),
        ]
        for cls, char in _hierarchy:
            try:
                if issubclass(dt, cls):
                    return dtype(char)
            except TypeError:
                pass

    return dtype('d')

class ndarray:
    __slots__ = ('_data', '_shape', '_dtype', '_order')

    # Class-level attributes required by numpy introspection
    __array_priority__    = 0.0
    __array_finalize__    = None
    __array_wrap__        = None
    __array_namespace__   = None
    __array_ufunc__       = None
    __array_function__    = None
    __copy__              = None
    __deepcopy__          = None
    byteswap              = None
    choose                = None
    compress              = None
    cumprod               = None
    cumsum                = None
    dump                  = None
    dumps                 = None
    fill                  = None
    getfield              = None
    setfield              = None
    setflags              = None
    item                  = None
    put                   = None
    repeat                = None
    resize                = None
    searchsorted          = None
    swapaxes              = None
    take                  = None
    to_device             = None
    argpartition          = None
    partition             = None

    def __init__(self, shape, dtype=None, buffer=None, offset=0,
                 strides=None, order='C'):
        if isinstance(shape, int):
            shape = (shape,)
        self._shape = tuple(shape)
        self._dtype = _coerce_dtype(dtype) or float64_dtype
        n = _prod(self._shape)
        if buffer is not None:
            flat = _flat_list(buffer)[:n]
            flat.extend([self._dtype.type() for _ in range(n - len(flat))])
            self._data = flat
        else:
            zero = 0 if self._dtype.kind not in ('U', 'S', 'O') else \
                   ('' if self._dtype.kind == 'U' else b'' if self._dtype.kind == 'S' else None)
            self._data = [zero] * n
        self._order = order

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, new_shape):
        if _prod(new_shape) != _prod(self._shape):
            raise ValueError('cannot reshape without changing size')
        self._shape = tuple(new_shape)

    @property
    def data(self): return memoryview(bytes(int(v) & 0xFF for v in self._data))
    @property
    def dtype(self): return self._dtype
    @property
    def ndim(self): return len(self._shape)
    @property
    def size(self): return _prod(self._shape)
    @property
    def itemsize(self): return self._dtype.itemsize
    @property
    def nbytes(self): return self.size * self.itemsize
    @property
    def strides(self):
        s = []
        acc = self.itemsize
        for d in reversed(self._shape):
            s.append(acc)
            acc *= d
        return tuple(reversed(s))
    @property
    def T(self):
        if self.ndim < 2:
            return self
        return self.transpose()
    @property
    def flat(self): return iter(self._data)
    @property
    def flags(self): return flagsobj()
    @property
    def real(self):
        if self._dtype.kind in ('F', 'D', 'G'):
            return _array_from_flat([x.real for x in self._data], self._shape, float64_dtype)
        return self.copy()
    @property
    def imag(self):
        if self._dtype.kind in ('F', 'D', 'G'):
            return _array_from_flat([x.imag for x in self._data], self._shape, float64_dtype)
        return zeros(self._shape, dtype=float64_dtype)

    def _multi_to_flat(self, idx):
        if not isinstance(idx, tuple):
            idx = (idx,)
        flat = 0
        stride = 1
        for d, i in zip(reversed(self._shape), reversed(idx)):
            if i < 0:
                i += d
            flat += i * stride
            stride *= d
        return flat

    def __getitem__(self, idx):
        if self.ndim == 0:
            return self._data[0]
        if isinstance(idx, tuple):
            expanded = []
            for i in idx:
                if i is Ellipsis:
                    n_missing = self.ndim - sum(1 for x in idx if x is not Ellipsis and x is not None)
                    expanded.extend([slice(None)] * n_missing)
                else:
                    expanded.append(i)
            idx = tuple(expanded)
        if idx is Ellipsis:
            idx = tuple(slice(None) for _ in range(self.ndim))
        if isinstance(idx, (int, slice)):
            idx = (idx,)
        if isinstance(idx, tuple):
            shape = list(self._shape)
            ndim = len(shape)
            full_idx = list(idx) + [slice(None)] * (ndim - len(idx))
            out_shape = []
            axis_selectors = []
            for ax, sel in enumerate(full_idx[:ndim]):
                dim = shape[ax]
                if isinstance(sel, int):
                    s = sel if sel >= 0 else sel + dim
                    axis_selectors.append(('int', s, dim))
                elif isinstance(sel, slice):
                    indices = list(range(*sel.indices(dim)))
                    axis_selectors.append(('slice', indices, dim))
                    out_shape.append(len(indices))
                elif sel is None:
                    axis_selectors.append(('new', None, 0))
                    out_shape.append(1)
                else:
                    idxs = _to_indices(sel, tuple(shape))
                    return _array_from_flat(
                        [self._data[j] for j in idxs],
                        (len(idxs),), self._dtype)
            src_strides = []
            s = 1
            for dim in reversed(shape):
                src_strides.insert(0, s)
                s *= dim
            out_data = []
            def _recurse(ax_sel_idx, src_offset):
                if ax_sel_idx == len(axis_selectors):
                    out_data.append(self._data[src_offset])
                    return
                kind_info = axis_selectors[ax_sel_idx]
                ax_stride = src_strides[ax_sel_idx]
                if kind_info[0] == 'int':
                    _recurse(ax_sel_idx + 1, src_offset + kind_info[1] * ax_stride)
                elif kind_info[0] == 'slice':
                    for j in kind_info[1]:
                        _recurse(ax_sel_idx + 1, src_offset + j * ax_stride)
                else:
                    _recurse(ax_sel_idx + 1, src_offset)
            _recurse(0, 0)
            if not out_shape:
                return out_data[0] if out_data else None
            return _array_from_flat(out_data, tuple(out_shape), self._dtype)
        return self._data[idx]

    def __setitem__(self, idx, value):
        if isinstance(idx, tuple):
            expanded = []
            for i in idx:
                if i is Ellipsis:
                    n_missing = self.ndim - sum(1 for x in idx if x is not Ellipsis and x is not None)
                    expanded.extend([slice(None)] * n_missing)
                else:
                    expanded.append(i)
            idx = tuple(expanded)
        if idx is Ellipsis:
            idx = tuple(slice(None) for _ in range(self.ndim))
        if isinstance(idx, (int, slice)):
            idx = (idx,)
        if isinstance(idx, tuple):
            shape = self._shape
            if all(isinstance(i, int) for i in idx):
                flat = 0
                for dim, i in zip(self._shape, idx):
                    if i < 0:
                        i += dim
                    flat = flat * dim + i
                if isinstance(value, ndarray):
                    value = value.tolist()
                if isinstance(value, list):
                    for k, v in enumerate(value):
                        self._data[flat + k] = v
                else:
                    self._data[flat] = value
                return
        if isinstance(value, ndarray):
            flat_vals = value._data
        elif isinstance(value, (list, tuple)):
            flat_vals = _flat_list(value)
        else:
            flat_vals = None
        if isinstance(idx, slice):
            indices = range(*idx.indices(self.size))
            if flat_vals is None:
                for i in indices:
                    self._data[i] = value
            else:
                for i, v in zip(indices, flat_vals):
                    self._data[i] = v
        else:
            self._data[idx] = value

    def reshape(self, *args, order='C'):
        if len(args) == 1 and isinstance(args[0], (tuple, list)):
            new_shape = tuple(args[0])
        else:
            new_shape = tuple(args)
        known = [x for x in new_shape if x != -1]
        if len(known) < len(new_shape):
            infer = self.size // (_prod(known) or 1)
            new_shape = tuple(infer if x == -1 else x for x in new_shape)
        if _prod(new_shape) != self.size:
            raise ValueError(f'cannot reshape size {self.size} into shape {new_shape}')
        result = ndarray.__new__(ndarray)
        result._data = list(self._data)
        result._shape = new_shape
        result._dtype = self._dtype
        result._order = order
        return result

    def flatten(self, order='C'):
        return _array_from_flat(list(self._data), (self.size,), self._dtype)

    def ravel(self, order='C'):
        return self.flatten(order)

    def transpose(self, axes=None):
        if self.ndim < 2:
            return self.copy()
        if axes is None:
            axes = tuple(reversed(range(self.ndim)))
        new_shape = tuple(self._shape[a] for a in axes)
        result_data = [None] * self.size
        old_strides = []
        s = 1
        for d in reversed(self._shape):
            old_strides.append(s)
            s *= d
        old_strides = list(reversed(old_strides))
        new_strides = []
        s = 1
        for d in reversed(new_shape):
            new_strides.append(s)
            s *= d
        new_strides = list(reversed(new_strides))
        for i in range(self.size):
            old_idx = []
            tmp = i
            for d in self._shape:
                old_idx.append(tmp % d)
                tmp //= d
            old_idx = list(reversed(old_idx))
            new_idx = [old_idx[a] for a in axes]
            new_flat = sum(ni * ns for ni, ns in zip(new_idx, new_strides))
            result_data[new_flat] = self._data[i]
        return _array_from_flat(result_data, new_shape, self._dtype)

    def copy(self, order='C'):
        return _array_from_flat(list(self._data), self._shape, self._dtype)

    def astype(self, dt, order='K', casting='unsafe', subok=True, copy=True):
        dt = _coerce_dtype(dt)
        conv = dt.type
        return _array_from_flat([conv(x) for x in self._data], self._shape, dt)

    def view(self, dtype=None):
        if dtype is None:
            return self.copy()
        return self.astype(dtype)

    def tolist(self):
        if self.ndim == 0:
            return self._data[0]
        def nest(data, shape):
            if len(shape) == 1:
                return list(data[:shape[0]])
            sub = _prod(shape[1:])
            return [nest(data[i * sub:(i + 1) * sub], shape[1:])
                    for i in range(shape[0])]
        return nest(self._data, self._shape)

    def tobytes(self, order='C'):
        fmt = self._dtype.str[1]
        try:
            return _struct.pack('=' + fmt * self.size, *self._data)
        except Exception:
            return bytes(self.size * self._dtype.itemsize)

    def tofile(self, fid, sep='', format='%s'):
        raise NotImplementedError('tofile() requires a filesystem')

    def fill(self, value):
        for i in range(self.size):
            self._data[i] = value

    def sum(self, axis=None, dtype=None, out=None, keepdims=False):
        if axis is None:
            return sum(self._data)
        return _reduce_axis(self, sum, axis, keepdims)

    def prod(self, axis=None, keepdims=False):
        if axis is None:
            return _prod(self._data)
        return _reduce_axis(self, _prod, axis, keepdims)

    def mean(self, axis=None, keepdims=False):
        if axis is None:
            return sum(self._data) / self.size
        return _reduce_axis(self, lambda d: sum(d) / len(d), axis, keepdims)

    def std(self, axis=None, ddof=0, keepdims=False):
        if axis is None:
            m = sum(self._data) / self.size
            return math.sqrt(sum((x - m) ** 2 for x in self._data) / (self.size - ddof))
        return _reduce_axis(self, lambda d: (lambda m: math.sqrt(sum((x-m)**2 for x in d)/(len(d)-ddof)))(sum(d)/len(d)), axis, keepdims)

    def var(self, axis=None, ddof=0, keepdims=False):
        if axis is None:
            m = sum(self._data) / self.size
            return sum((x - m) ** 2 for x in self._data) / (self.size - ddof)
        return _reduce_axis(self, lambda d: (lambda m: sum((x-m)**2 for x in d)/(len(d)-ddof))(sum(d)/len(d)), axis, keepdims)

    def min(self, axis=None, keepdims=False):
        if axis is None:
            return min(self._data)
        return _reduce_axis(self, min, axis, keepdims)

    def max(self, axis=None, keepdims=False):
        if axis is None:
            return max(self._data)
        return _reduce_axis(self, max, axis, keepdims)

    def argmin(self, axis=None):
        if axis is None:
            return self._data.index(min(self._data))
        raise NotImplementedError

    def argmax(self, axis=None):
        if axis is None:
            return self._data.index(max(self._data))
        raise NotImplementedError

    def argsort(self, axis=-1, kind=None, order=None):
        if axis is None or self.ndim <= 1:
            indexed = sorted(enumerate(self._data), key=lambda x: x[1])
            return _array_from_flat([i for i, _ in indexed], self._shape, int64_dtype)
        raise NotImplementedError('argsort with axis on multi-dim array')

    def sort(self, axis=-1, kind=None, order=None):
        if axis is None or self.ndim <= 1:
            self._data.sort()
        else:
            raise NotImplementedError

    def nonzero(self):
        indices = [i for i, v in enumerate(self._data) if v]
        if self.ndim == 1:
            return (_array_from_flat(indices, (len(indices),), int64_dtype),)
        result = [[] for _ in range(self.ndim)]
        for flat_i in indices:
            multi = _unravel(flat_i, self._shape)
            for dim, idx in enumerate(multi):
                result[dim].append(idx)
        return tuple(_array_from_flat(r, (len(r),), int64_dtype) for r in result)

    def diagonal(self, offset=0, axis1=0, axis2=1):
        if self.ndim != 2:
            raise NotImplementedError
        rows, cols = self._shape
        diag = []
        for i in range(min(rows, cols - offset) if offset >= 0 else min(rows + offset, cols)):
            r = i if offset >= 0 else i - offset
            c = i + offset if offset >= 0 else i
            diag.append(self._data[r * cols + c])
        return _array_from_flat(diag, (len(diag),), self._dtype)

    def trace(self, offset=0):
        return sum(self.diagonal(offset)._data)

    def dot(self, b):
        return dot(self, b)

    def __matmul__(self, other):
        return matmul(self, other)

    def __iter__(self):
        if self.ndim == 0:
            raise TypeError('iteration over a 0-d array')
        sub_size = _prod(self._shape[1:]) if self.ndim > 1 else 1
        sub_shape = self._shape[1:] if self.ndim > 1 else ()
        for i in range(self._shape[0]):
            chunk = self._data[i * sub_size:(i + 1) * sub_size]
            if sub_shape:
                yield _array_from_flat(chunk, sub_shape, self._dtype)
            else:
                yield chunk[0]

    def __len__(self):
        if self.ndim == 0:
            raise TypeError('len() of unsized object')
        return self._shape[0]

    def __repr__(self):
        return f'array({self.tolist()!r}, dtype={self._dtype.name})'

    def __str__(self):
        return str(self.tolist())

    def __array__(self, dt=None):
        if dt is None:
            return self
        return self.astype(dt)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        return NotImplemented

    def __array_function__(self, func, types, args, kwargs):
        return NotImplemented

    def __add__(self, o): return _ew(operator.add, self, o)
    def __radd__(self, o): return _ew(operator.add, o, self)
    def __sub__(self, o): return _ew(operator.sub, self, o)
    def __rsub__(self, o): return _ew(operator.sub, o, self)
    def __mul__(self, o): return _ew(operator.mul, self, o)
    def __rmul__(self, o): return _ew(operator.mul, o, self)
    def __truediv__(self, o): return _ew(operator.truediv, self, o)
    def __rtruediv__(self, o): return _ew(operator.truediv, o, self)
    def __floordiv__(self, o): return _ew(operator.floordiv, self, o)
    def __mod__(self, o): return _ew(operator.mod, self, o)
    def __pow__(self, o): return _ew(operator.pow, self, o)
    def __neg__(self): return _array_from_flat([-x for x in self._data], self._shape, self._dtype)
    def __pos__(self): return self.copy()
    def __abs__(self): return _array_from_flat([abs(x) for x in self._data], self._shape, self._dtype)
    def __invert__(self): return _array_from_flat([~x for x in self._data], self._shape, self._dtype)
    def __and__(self, o): return _ew(operator.and_, self, o)
    def __or__(self, o): return _ew(operator.or_, self, o)
    def __xor__(self, o): return _ew(operator.xor, self, o)
    def __lshift__(self, o): return _ew(operator.lshift, self, o)
    def __rshift__(self, o): return _ew(operator.rshift, self, o)
    def __lt__(self, o): return _ew(operator.lt, self, o)
    def __le__(self, o): return _ew(operator.le, self, o)
    def __gt__(self, o): return _ew(operator.gt, self, o)
    def __ge__(self, o): return _ew(operator.ge, self, o)
    def __eq__(self, o): return _ew(operator.eq, self, o)
    def __ne__(self, o): return _ew(operator.ne, self, o)
    def __bool__(self):
        if self.size == 1:
            return bool(self._data[0])
        raise ValueError('The truth value of an array with more than one element is ambiguous.')
    def __float__(self):
        if self.size == 1:
            return float(self._data[0])
        raise TypeError
    def __int__(self):
        if self.size == 1:
            return int(self._data[0])
        raise TypeError
    def __complex__(self):
        if self.size == 1:
            return complex(self._data[0])
        raise TypeError

    def all(self, axis=None, keepdims=False):
        if axis is None:
            return all(self._data)
        return _reduce_axis(self, all, axis, keepdims)

    def any(self, axis=None, keepdims=False):
        if axis is None:
            return any(self._data)
        return _reduce_axis(self, any, axis, keepdims)

    def clip(self, a_min=None, a_max=None, out=None):
        def _clp(x):
            if a_min is not None and x < a_min:
                return a_min
            if a_max is not None and x > a_max:
                return a_max
            return x
        return _array_from_flat([_clp(x) for x in self._data], self._shape, self._dtype)

    def round(self, decimals=0):
        return _array_from_flat([round(x, decimals) for x in self._data], self._shape, self._dtype)

    def squeeze(self, axis=None):
        new_shape = tuple(d for d in self._shape if d != 1)
        return _array_from_flat(list(self._data), new_shape or (1,), self._dtype)

    def conj(self):
        if self._dtype.kind in ('F', 'D', 'G'):
            return _array_from_flat([x.conjugate() for x in self._data], self._shape, self._dtype)
        return self.copy()

    conjugate = conj

    def __contains__(self, item):
        return item in self._data

    def __hash__(self):
        return id(self)

    __array_priority__ = 0.0
    __array_finalize__ = None
    __array_wrap__ = None

    @property
    def __array_interface__(self):
        return {'shape': self._shape, 'typestr': self._dtype.str,
                'data': (id(self), False), 'version': 3}
    @property
    def __array_struct__(self): return None
    @property
    def mT(self): return self.T
    @property
    def base(self): return None
    @property
    def ctypes(self): return None

    @classmethod
    def __class_getitem__(cls, item):
        return types.GenericAlias(cls, (item,) if not isinstance(item, tuple) else item)

    def __dlpack__(self, stream=None):
        raise NotImplementedError("__dlpack__ not supported in pure-Python shim")
    def __dlpack_device__(self): return (1, 0)

    def __reduce__(self):
        return (_reconstruct, (ndarray, (0,), b'b'), self.__getstate__())
    def __reduce_ex__(self, protocol): return self.__reduce__()
    def __getstate__(self):
        return (1, self._shape, self._dtype, False,
                bytes(int(v) & 0xFF for v in self._data))
    def __setstate__(self, state):
        _version, shape, dt, _fortran, data = state
        self._shape = tuple(shape)
        self._dtype = dtype(dt) if not isinstance(dt, dtype) else dt
        self._data = list(data)

    def argpartition(self, kth, axis=-1, kind='introselect', order=None):
        import builtins
        flat = list(self._data)
        idx = builtins.sorted(range(len(flat)), key=lambda i: flat[i])
        return ndarray.__new__(ndarray)
    def partition(self, kth, axis=-1, kind='introselect', order=None):
        self._data.sort()

    # Instance method implementations for class-level stubs
    def __array_namespace__(self, *, api_version=None): return None
    def __copy__(self): return self.copy()
    def __deepcopy__(self, memo=None): return self.copy()
    def byteswap(self, inplace=False): return self if inplace else self.copy()
    def choose(self, choices, out=None, mode='raise'):
        import builtins
        result = [choices[int(v)] for v in self._data]
        flat = []
        for c in result:
            a = array(c)
            flat.extend(a._data)
        return _array_from_flat(flat, self._shape, self._dtype)
    def compress(self, condition, axis=None, out=None):
        cond = [bool(v) for v in array(condition)._data]
        flat = [v for v, c in zip(self._data, cond) if c]
        return _array_from_flat(flat, (len(flat),), self._dtype)
    def cumprod(self, axis=None, dtype=None, out=None):
        import functools, operator
        flat = list(self._data)
        result = []
        acc = 1
        for v in flat:
            acc *= v
            result.append(acc)
        return _array_from_flat(result, (len(result),), self._dtype)
    def cumsum(self, axis=None, dtype=None, out=None):
        flat = list(self._data)
        result = []
        acc = 0
        for v in flat:
            acc += v
            result.append(acc)
        return _array_from_flat(result, (len(result),), self._dtype)
    def dump(self, file):
        import pickle
        with open(file, 'wb') as f: pickle.dump(self, f)
    def dumps(self):
        import pickle
        return pickle.dumps(self)
    def getfield(self, dtype, offset=0): return self.view(dtype)
    def setfield(self, val, dtype, offset=0): pass
    def setflags(self, write=None, align=None, uic=None): pass
    def item(self, *args):
        if not args:
            if self.size != 1:
                raise ValueError("can only convert an array of size 1 to a Python scalar")
            return self._data[0]
        return self._data[args[0]]
    def put(self, indices, values, mode='raise'):
        idx_list = array(indices)._data
        val_list = array(values)._data
        for i, v in zip(idx_list, val_list):
            self._data[int(i)] = v
    def repeat(self, repeats, axis=None):
        flat = []
        for v in self._data:
            flat.extend([v] * (repeats if isinstance(repeats, int) else repeats[0]))
        return _array_from_flat(flat, (len(flat),), self._dtype)
    def resize(self, new_shape, refcheck=True):
        if isinstance(new_shape, int):
            new_shape = (new_shape,)
        new_size = _prod(new_shape)
        old = list(self._data)
        new_data = (old * (new_size // len(old) + 1))[:new_size] if old else [0] * new_size
        self._data = new_data
        self._shape = tuple(new_shape)
    def searchsorted(self, v, side='left', sorter=None):
        import bisect
        flat = sorted(self._data) if sorter is None else [self._data[i] for i in sorter]
        fn = bisect.bisect_left if side == 'left' else bisect.bisect_right
        vv = array(v)
        result = [fn(flat, x) for x in vv._data]
        return _array_from_flat(result, vv._shape, dtype('int64'))
    def swapaxes(self, axis1, axis2):
        if self.ndim < 2:
            return self.copy()
        new_shape = list(self._shape)
        new_shape[axis1], new_shape[axis2] = new_shape[axis2], new_shape[axis1]
        return _array_from_flat(list(self._data), tuple(new_shape), self._dtype)
    def take(self, indices, axis=None, out=None, mode='raise'):
        idx = array(indices)
        result = [self._data[int(i)] for i in idx._data]
        return _array_from_flat(result, idx._shape, self._dtype)
    def to_device(self, device, *, stream=None): return self

def _to_indices(idx_arr, shape):
    if isinstance(idx_arr, ndarray):
        if idx_arr._dtype.kind == 'b':
            return [i for i, v in enumerate(idx_arr._data) if v]
        return [int(x) for x in idx_arr._data]
    return list(idx_arr)

def _unwrap(v):
    if isinstance(v, generic):
        return v._value
    return v

def _array_from_flat(data, shape, dt):
    arr = ndarray.__new__(ndarray)
    arr._data = [_unwrap(v) for v in data]
    arr._shape = tuple(shape)
    arr._dtype = dt
    arr._order = 'C'
    return arr

def _unravel(flat, shape):
    idx = []
    for d in reversed(shape):
        idx.append(flat % d)
        flat //= d
    return tuple(reversed(idx))

def _broadcast_shapes(s1, s2):
    r = max(len(s1), len(s2))
    s1 = (1,) * (r - len(s1)) + s1
    s2 = (1,) * (r - len(s2)) + s2
    out = []
    for a, b in zip(s1, s2):
        if a == b: out.append(a)
        elif a == 1: out.append(b)
        elif b == 1: out.append(a)
        else: raise ValueError(f'cannot broadcast shapes {s1} and {s2}')
    return tuple(out)

def _broadcast_index(flat, src_shape, dst_shape):
    r = len(dst_shape)
    ps = (1,) * (r - len(src_shape)) + src_shape
    dst_idx = _unravel(flat, dst_shape)
    src_idx = tuple(0 if ps[i] == 1 else dst_idx[i] for i in range(r))
    src_flat = 0
    for i, d in zip(src_idx, ps):
        src_flat = src_flat * d + i
    return src_flat

def _ew(op, a, b):
    if not isinstance(a, ndarray):
        a_data, a_shape = [a], ()
    else:
        a_data, a_shape = a._data, a._shape
    if not isinstance(b, ndarray):
        b_data, b_shape = [b], ()
    else:
        b_data, b_shape = b._data, b._shape
    if not a_shape and not b_shape:
        return op(a_data[0], b_data[0])
    out_shape = _broadcast_shapes(a_shape or (1,), b_shape or (1,))
    n = _prod(out_shape)
    out = []
    for i in range(n):
        ai = _broadcast_index(i, a_shape or (1,), out_shape) if a_shape else 0
        bi = _broadcast_index(i, b_shape or (1,), out_shape) if b_shape else 0
        av = a_data[ai] if a_shape else a_data[0]
        bv = b_data[bi] if b_shape else b_data[0]
        try:
            out.append(op(av, bv))
        except (ZeroDivisionError, ArithmeticError):
            out.append(float('nan') if isinstance(av, float) else 0)
    dt = (a if isinstance(a, ndarray) else b)._dtype
    return _array_from_flat(out, out_shape, dt)

def _reduce_axis(arr, fn, axis, keepdims):
    if axis < 0:
        axis += arr.ndim
    sub = _prod(arr._shape[axis + 1:]) if axis + 1 < arr.ndim else 1
    stride = sub * arr._shape[axis]
    new_shape = arr._shape[:axis] + (arr._shape[axis + 1:] if not keepdims else (1,) + arr._shape[axis + 1:])
    outer = _prod(arr._shape[:axis])
    out = []
    for o in range(outer):
        for s in range(sub):
            chunk = [arr._data[o * stride + j * sub + s] for j in range(arr._shape[axis])]
            out.append(fn(chunk))
    return _array_from_flat(out, new_shape or (1,), arr._dtype)

class flagsobj:
    def __init__(self):
        self.c_contiguous = True
        self.f_contiguous = False
        self.owndata = True
        self.writeable = True
        self.aligned = True
        self.writebackifcopy = False
        self.updateifcopy = False
    def __getitem__(self, k):
        return getattr(self, k.lower().replace(' ', '_').replace('-', '_'), False)
    def __setitem__(self, k, v):
        setattr(self, k.lower().replace(' ', '_').replace('-', '_'), v)
    def __repr__(self):
        return (f'  C_CONTIGUOUS : {self.c_contiguous}\n'
                f'  F_CONTIGUOUS : {self.f_contiguous}\n'
                f'  OWNDATA : {self.owndata}\n'
                f'  WRITEABLE : {self.writeable}\n'
                f'  ALIGNED : {self.aligned}')

# ===========================================================================
# Section 4 – Array creation functions
# ===========================================================================
def array(obj, dtype=None, *, copy=True, order='K', subok=False, ndmin=0, like=None):
    if isinstance(obj, ndarray) and not copy and dtype is None:
        return obj
    if isinstance(obj, ndarray):
        data = list(obj._data)
        shape = obj._shape
    elif isinstance(obj, (list, tuple)):
        shape = _nested_shape(obj)
        data = _flat_list(obj)
    elif hasattr(obj, '__iter__'):
        data = list(obj)
        shape = (len(data),)
    else:
        data = [obj]
        shape = ()
    while len(shape) < ndmin:
        shape = (1,) + shape
    dt = _coerce_dtype(dtype)
    if dt is None:
        if data:
            sample = data[0]
            if isinstance(sample, bool):
                dt = bool_dtype
            elif isinstance(sample, int):
                dt = int64_dtype
            elif isinstance(sample, float):
                dt = float64_dtype
            elif isinstance(sample, complex):
                dt = complex128_dtype
            elif isinstance(sample, bytes):
                dt = dtype('S')
            elif isinstance(sample, str):
                dt = dtype('U')
            else:
                dt = object_dtype
        else:
            dt = float64_dtype
    conv = dt.type
    try:
        data = [conv(x) for x in data]
    except (TypeError, ValueError):
        pass
    return _array_from_flat(data, shape, dt)

def _py_zero(dt):
    k = dt.kind
    if k == 'b': return False
    if k in ('i','u'): return 0
    if k in ('f',): return 0.0
    if k == 'c': return 0j
    if k == 'S': return b''
    if k == 'U': return ''
    return None

def _py_one(dt):
    k = dt.kind
    if k == 'b': return True
    if k in ('i','u'): return 1
    if k in ('f',): return 1.0
    if k == 'c': return 1+0j
    return 1

def zeros(shape, dtype=None, order='C', *, like=None, device=None):
    if isinstance(shape, int):
        shape = (shape,)
    shape = tuple(shape)
    dt = _coerce_dtype(dtype) or float64_dtype
    n = _prod(shape)
    return _array_from_flat([_py_zero(dt)] * n, shape, dt)

def ones(shape, dtype=None, order='C', *, like=None, device=None):
    dt = _coerce_dtype(dtype) or float64_dtype
    if isinstance(shape, int):
        shape = (shape,)
    shape = tuple(shape)
    return _array_from_flat([_py_one(dt)] * _prod(shape), shape, dt)

def empty(shape, dtype=None, order='C', *, like=None, device=None):
    return zeros(shape, dtype=dtype, order=order)

def full(shape, fill_value, dtype=None, order='C', *, like=None):
    dt = _coerce_dtype(dtype)
    if isinstance(shape, int):
        shape = (shape,)
    shape = tuple(shape)
    if dt is None:
        if isinstance(fill_value, bool):
            dt = bool_dtype
        elif isinstance(fill_value, int):
            dt = int64_dtype
        elif isinstance(fill_value, float):
            dt = float64_dtype
        elif isinstance(fill_value, complex):
            dt = complex128_dtype
        else:
            dt = object_dtype
    return _array_from_flat([dt.type(fill_value)] * _prod(shape), shape, dt)

def zeros_like(a, dtype=None, order='K', subok=True, shape=None):
    dt = _coerce_dtype(dtype) or (a._dtype if isinstance(a, ndarray) else float64_dtype)
    sh = shape or (a._shape if isinstance(a, ndarray) else ())
    return zeros(sh, dtype=dt)

def ones_like(a, dtype=None, order='K', subok=True, shape=None):
    dt = _coerce_dtype(dtype) or (a._dtype if isinstance(a, ndarray) else float64_dtype)
    sh = shape or (a._shape if isinstance(a, ndarray) else ())
    return ones(sh, dtype=dt)

def empty_like(a, dtype=None, order='K', subok=True, shape=None):
    return zeros_like(a, dtype=dtype, shape=shape)

def full_like(a, fill_value, dtype=None, order='K', subok=True, shape=None):
    dt = _coerce_dtype(dtype) or (a._dtype if isinstance(a, ndarray) else float64_dtype)
    sh = shape or (a._shape if isinstance(a, ndarray) else ())
    return full(sh, fill_value, dtype=dt)

def arange(start, stop=None, step=1, dtype=None, *, like=None, device=None):
    if stop is None:
        start, stop = 0, start
    data = []
    x = start
    if step > 0:
        while x < stop:
            data.append(x)
            x += step
    else:
        while x > stop:
            data.append(x)
            x += step
    dt = _coerce_dtype(dtype)
    if dt is None:
        if isinstance(step, float) or isinstance(start, float) or isinstance(stop, float):
            dt = float64_dtype
        else:
            dt = int64_dtype
    return _array_from_flat([dt.type(x) for x in data], (len(data),), dt)

def linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None, axis=0):
    n = int(num)
    if n <= 0:
        data = []
        step = float('nan')
    elif n == 1:
        data = [float(start)]
        step = 0.0
    else:
        step = (stop - start) / (n - 1 if endpoint else n)
        data = [start + i * step for i in range(n)]
        if endpoint and n > 1:
            data[-1] = stop
    dt = _coerce_dtype(dtype) or float64_dtype
    arr = _array_from_flat([dt.type(x) for x in data], (len(data),), dt)
    return (arr, step) if retstep else arr

def logspace(start, stop, num=50, endpoint=True, base=10.0, dtype=None, axis=0):
    lin = linspace(start, stop, num=num, endpoint=endpoint)
    data = [base ** x for x in lin._data]
    dt = _coerce_dtype(dtype) or float64_dtype
    return _array_from_flat([dt.type(x) for x in data], (len(data),), dt)

def eye(N, M=None, k=0, dtype=None, order='C', *, like=None):
    M = M or N
    dt = _coerce_dtype(dtype) or float64_dtype
    data = [dt.type(1 if j - i == k else 0) for i in range(N) for j in range(M)]
    return _array_from_flat(data, (N, M), dt)

def identity(n, dtype=None, *, like=None):
    return eye(n, dtype=dtype)

def diag(v, k=0):
    if isinstance(v, ndarray) and v.ndim == 1:
        n = v.size + abs(k)
        dt = v._dtype
        data = [dt.type(0)] * (n * n)
        for i in range(v.size):
            r = i if k >= 0 else i - k
            c = i + k if k >= 0 else i
            data[r * n + c] = v._data[i]
        return _array_from_flat(data, (n, n), dt)
    elif isinstance(v, ndarray) and v.ndim == 2:
        return v.diagonal(k)
    raise TypeError

def diagflat(v, k=0):
    return diag(array(v).flatten(), k)

def tri(N, M=None, k=0, dtype=None, *, like=None):
    M = M or N
    dt = _coerce_dtype(dtype) or float64_dtype
    data = [dt.type(1 if j <= i + k else 0) for i in range(N) for j in range(M)]
    return _array_from_flat(data, (N, M), dt)

def tril(m, k=0):
    a = array(m)
    rows, cols = a._shape[-2], a._shape[-1]
    data = list(a._data)
    for i in range(rows):
        for j in range(cols):
            if j > i + k:
                data[i * cols + j] = a._dtype.type(0)
    return _array_from_flat(data, a._shape, a._dtype)

def triu(m, k=0):
    a = array(m)
    rows, cols = a._shape[-2], a._shape[-1]
    data = list(a._data)
    for i in range(rows):
        for j in range(cols):
            if j < i + k:
                data[i * cols + j] = a._dtype.type(0)
    return _array_from_flat(data, a._shape, a._dtype)

def concatenate(arrays, axis=0, out=None, *, dtype=None, casting='same_kind'):
    arrays = [array(a) for a in arrays]
    if not arrays:
        raise ValueError('need at least one array to concatenate')
    if axis < 0:
        axis += arrays[0].ndim
    dt = _coerce_dtype(dtype) or arrays[0]._dtype
    if arrays[0].ndim <= 1:
        data = []
        for a in arrays:
            data.extend(a._data)
        return _array_from_flat(data, (len(data),), dt)
    if axis == 0:
        data = []
        for a in arrays:
            data.extend(a._data)
        new_shape = (sum(a._shape[0] for a in arrays),) + arrays[0]._shape[1:]
        return _array_from_flat(data, new_shape, dt)
    if axis == 1:
        rows = arrays[0]._shape[0]
        cols_per = [a._shape[1] for a in arrays]
        total_cols = sum(cols_per)
        data = []
        for r in range(rows):
            for a in arrays:
                c = a._shape[1]
                data.extend(a._data[r * c:(r + 1) * c])
        return _array_from_flat(data, (rows, total_cols), dt)
    raise NotImplementedError(f'concatenate along axis={axis} not supported')

def stack(arrays, axis=0, out=None, *, dtype=None, casting='same_kind'):
    arrays = [array(a) for a in arrays]
    n = len(arrays)
    sh = arrays[0]._shape
    new_shape = sh[:axis] + (n,) + sh[axis:]
    expanded = [a.reshape(sh[:axis] + (1,) + sh[axis:]) for a in arrays]
    return concatenate(expanded, axis=axis)

def vstack(tup):
    return concatenate([atleast_2d(a) for a in tup], axis=0)

def hstack(tup):
    arrays = [atleast_1d(a) for a in tup]
    if arrays[0].ndim == 1:
        return concatenate(arrays, axis=0)
    return concatenate(arrays, axis=1)

def dstack(tup):
    return concatenate([atleast_3d(a) for a in tup], axis=2)

def atleast_1d(*arys):
    result = []
    for a in arys:
        a = array(a)
        if a.ndim == 0:
            a = a.reshape(1)
        result.append(a)
    return result[0] if len(result) == 1 else result

def atleast_2d(*arys):
    result = []
    for a in arys:
        a = array(a)
        if a.ndim == 0:
            a = a.reshape(1, 1)
        elif a.ndim == 1:
            a = a.reshape(1, a.size)
        result.append(a)
    return result[0] if len(result) == 1 else result

def atleast_3d(*arys):
    result = []
    for a in arys:
        a = array(a)
        if a.ndim == 0:
            a = a.reshape(1, 1, 1)
        elif a.ndim == 1:
            a = a.reshape(1, a.size, 1)
        elif a.ndim == 2:
            a = a.reshape(a._shape + (1,))
        result.append(a)
    return result[0] if len(result) == 1 else result

def split(ary, indices_or_sections, axis=0):
    ary = array(ary)
    n = ary._shape[axis]
    if isinstance(indices_or_sections, int):
        k = indices_or_sections
        indices = list(range(0, n, n // k))
    else:
        indices = [0] + list(indices_or_sections) + [n]
    result = []
    for i in range(len(indices) - 1):
        sl = [slice(None)] * ary.ndim
        sl[axis] = slice(indices[i], indices[i + 1])
        result.append(ary[tuple(sl)])
    return result

def hsplit(ary, indices_or_sections):
    return split(ary, indices_or_sections, axis=1 if array(ary).ndim > 1 else 0)

def vsplit(ary, indices_or_sections):
    return split(ary, indices_or_sections, axis=0)

def repeat(a, repeats, axis=None):
    a = array(a)
    if axis is None:
        data = []
        for x in a._data:
            n = repeats if isinstance(repeats, int) else repeats
            for _ in range(n):
                data.append(x)
        return _array_from_flat(data, (len(data),), a._dtype)
    raise NotImplementedError

def tile(A, reps):
    A = array(A)
    if isinstance(reps, int):
        reps = (reps,)
    result = A
    for axis, r in enumerate(reversed(reps)):
        result = concatenate([result] * r, axis=result.ndim - 1 - axis)
    return result

def where(condition, x=None, y=None):
    condition = array(condition)
    if x is None and y is None:
        return condition.nonzero()
    x = array(x) if not isinstance(x, ndarray) else x
    y = array(y) if not isinstance(y, ndarray) else y
    out_shape = condition._shape
    data = []
    for i, c in enumerate(condition._data):
        ci = _broadcast_index(i, condition._shape, out_shape) if condition._shape else 0
        xi = _broadcast_index(i, x._shape, out_shape) if x._shape else 0
        yi = _broadcast_index(i, y._shape, out_shape) if y._shape else 0
        data.append(x._data[xi] if condition._data[ci] else y._data[yi])
    return _array_from_flat(data, out_shape, x._dtype)

def nonzero(a):
    return array(a).nonzero()

def argwhere(a):
    a = array(a)
    result = [list(_unravel(i, a._shape)) for i, v in enumerate(a._data) if v]
    if not result:
        return zeros((0, a.ndim), dtype=int64_dtype)
    return _array_from_flat([x for row in result for x in row], (len(result), a.ndim), int64_dtype)

def ravel(a, order='C'):
    return array(a).flatten()

def take(a, indices, axis=None, out=None, mode='raise'):
    a = array(a)
    indices = array(indices)
    data = [a._data[int(i)] for i in indices._data]
    return _array_from_flat(data, indices._shape, a._dtype)

def put(a, ind, v, mode='raise'):
    a = array(a)
    for i, vi in zip(_flat_list(ind), _flat_list(v)):
        a._data[int(i)] = vi

def copyto(dst, src, casting='same_kind', where=True):
    src_arr = array(src)
    n = dst.size
    sd = src_arr._data
    if len(sd) == 0:
        return
    if len(sd) == 1:
        dst._data[:] = [sd[0]] * n
    elif len(sd) >= n:
        dst._data[:] = list(sd[:n])
    else:
        dst._data[:] = [sd[i % len(sd)] for i in range(n)]

def asarray(a, dtype=None, order=None, *, like=None):
    return array(a, dtype=dtype, copy=False)

def asanyarray(a, dtype=None, order=None, like=None):
    return asarray(a, dtype=dtype)

def ascontiguousarray(a, dtype=None, *, like=None):
    return array(a, dtype=dtype, order='C')

def asfortranarray(a, dtype=None, *, like=None):
    return array(a, dtype=dtype, order='F')

def frombuffer(buffer, dtype='float64', count=-1, offset=0):
    dt = _coerce_dtype(dtype)
    raw = bytes(buffer)[offset:]
    fmt = dt.str[1]
    sz = dt.itemsize
    n = len(raw) // sz if count < 0 else count
    data = []
    for i in range(n):
        chunk = raw[i * sz:(i + 1) * sz]
        try:
            val = _struct.unpack('=' + fmt, chunk)[0]
        except Exception:
            val = 0
        data.append(dt.type(val))
    return _array_from_flat(data, (n,), dt)

def fromiter(iterable, dtype, count=-1):
    dt = _coerce_dtype(dtype)
    data = list(iterable) if count < 0 else list(itertools.islice(iterable, count))
    return _array_from_flat([dt.type(x) for x in data], (len(data),), dt)

def fromstring(string, dtype='float64', count=-1, *, sep=''):
    dt = _coerce_dtype(dtype)
    if sep:
        data = [dt.type(x) for x in string.split(sep)]
    else:
        return frombuffer(string.encode() if isinstance(string, str) else string, dtype=dt, count=count)
    return _array_from_flat(data, (len(data),), dt)

def frompyfunc(func, nin, nout, *, identity=None):
    class _pyfunc_ufunc:
        def __call__(self, *args, **kw):
            arrays = [array(a) for a in args]
            out_shape = functools.reduce(_broadcast_shapes,
                                         [a._shape for a in arrays],
                                         ())
            n = _prod(out_shape) if out_shape else 1
            results = [[] for _ in range(nout)]
            for i in range(n):
                in_vals = []
                for a in arrays:
                    idx = _broadcast_index(i, a._shape, out_shape) if a._shape else 0
                    in_vals.append(a._data[idx] if a._shape else a._data[0])
                out_vals = func(*in_vals)
                if nout == 1:
                    out_vals = (out_vals,)
                for k, v in enumerate(out_vals):
                    results[k].append(v)
            if nout == 1:
                return _array_from_flat(results[0], out_shape or (1,), object_dtype)
            return tuple(_array_from_flat(r, out_shape or (1,), object_dtype) for r in results)
        @property
        def nin(self): return nin
        @property
        def nout(self): return nout
        @property
        def nargs(self): return nin + nout
    return _pyfunc_ufunc()

# ===========================================================================
# Section 5 – Matrix / linear algebra helpers
# ===========================================================================
def dot(a, b, out=None):
    a, b = array(a), array(b)
    if a.ndim == 1 and b.ndim == 1:
        result = sum(x * y for x, y in zip(a._data, b._data))
        out_dt = a._dtype if a._dtype.kind == 'f' or a._dtype.kind == 'c' else b._dtype
        return _array_from_flat([result], (), out_dt)
    if a.ndim == 2 and b.ndim == 2:
        rows, inner = a._shape
        inner2, cols = b._shape
        if inner != inner2:
            raise ValueError(f'shapes {a._shape} and {b._shape} not aligned')
        data = []
        for i in range(rows):
            for j in range(cols):
                data.append(sum(a._data[i * inner + k] * b._data[k * cols + j]
                                for k in range(inner)))
        result = _array_from_flat(data, (rows, cols), a._dtype)
        if out is not None:
            out._data[:] = result._data
        return result
    if a.ndim == 2 and b.ndim == 1:
        rows, cols = a._shape
        data = [sum(a._data[i * cols + k] * b._data[k] for k in range(cols))
                for i in range(rows)]
        return _array_from_flat(data, (rows,), a._dtype)
    if a.ndim == 0 or b.ndim == 0:
        return _ew(operator.mul, a, b)
    raise NotImplementedError(f'dot for ndim {a.ndim},{b.ndim}')

def matmul(a, b, out=None, **kw):
    return dot(a, b, out)

def vdot(a, b):
    a, b = array(a).flatten(), array(b).flatten()
    return sum(x.conjugate() * y if isinstance(x, complex) else x * y
               for x, y in zip(a._data, b._data))

def inner(a, b):
    a, b = array(a), array(b)
    if a.ndim == 1 and b.ndim == 1:
        return sum(x * y for x, y in zip(a._data, b._data))
    raise NotImplementedError('inner for ndim > 1')

def outer(a, b, out=None):
    a, b = array(a).flatten(), array(b).flatten()
    data = [x * y for x in a._data for y in b._data]
    return _array_from_flat(data, (a.size, b.size), a._dtype)

def tensordot(a, b, axes=2):
    a, b = array(a), array(b)
    if isinstance(axes, int):
        axes_a = list(range(-axes, 0))
        axes_b = list(range(0, axes))
    else:
        axes_a, axes_b = axes
    if len(axes_a) == a.ndim and len(axes_b) == b.ndim:
        return sum(x * y for x, y in zip(a._data, b._data))
    na = _prod([a._shape[i] for i in axes_a])
    nb = _prod([b._shape[i] for i in axes_b])
    a2 = a.reshape(-1, na)
    b2 = b.reshape(nb, -1)
    return dot(a2, b2)

def kron(a, b):
    a, b = array(a), array(b)
    data = [x * y for x in a._data for y in b._data]
    shape = tuple(sa * sb for sa, sb in zip(a._shape, b._shape))
    return _array_from_flat(data, shape, a._dtype)

def cross(a, b, axisa=-1, axisb=-1, axisc=-1, axis=None):
    a, b = array(a).flatten(), array(b).flatten()
    if a.size == 3 and b.size == 3:
        ax, ay, az = a._data
        bx, by, bz = b._data
        return _array_from_flat([ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx],
                                (3,), float64_dtype)
    if a.size == 2 and b.size == 2:
        return a._data[0]*b._data[1] - a._data[1]*b._data[0]
    raise NotImplementedError

def einsum(subscripts, *operands, **kw):
    operands = [array(o) for o in operands]
    subscripts = subscripts.replace(' ', '')
    if '->' in subscripts:
        lhs, rhs = subscripts.split('->')
        inputs = lhs.split(',')
    else:
        inputs = [subscripts]
        rhs = None
    if len(inputs) == 1 and rhs is not None:
        a = operands[0]
        if a.ndim == 2 and set(inputs[0]) == set(rhs):
            if inputs[0] == rhs[::-1]:
                return a.transpose()
        return a
    if len(inputs) == 2 and rhs is not None:
        a, b = operands
        return dot(a, b)
    raise NotImplementedError(f'einsum pattern {subscripts!r} not supported')

def c_einsum(subscripts, *operands, **kw):
    return einsum(subscripts, *operands, **kw)

def trace(a, offset=0, axis1=0, axis2=1, dtype=None, out=None):
    return array(a).trace(offset)

def norm(x, ord=None, axis=None, keepdims=False):
    x = array(x)
    if axis is None:
        data = x._data
        if ord is None or ord == 2:
            return math.sqrt(sum(abs(v) ** 2 for v in data))
        if ord == 1:
            return sum(abs(v) for v in data)
        if ord == math.inf:
            return max(abs(v) for v in data)
        return sum(abs(v) ** ord for v in data) ** (1.0 / ord)
    raise NotImplementedError

def count_nonzero(a, axis=None, *, keepdims=False):
    a = array(a)
    if axis is None:
        return sum(1 for x in a._data if x)
    return _reduce_axis(a, lambda d: sum(1 for x in d if x), axis, keepdims)

def bincount(x, weights=None, minlength=0):
    x = array(x)
    max_val = max(int(v) for v in x._data) if x._data else -1
    n = max(max_val + 1, minlength)
    counts = [0] * n
    for v in x._data:
        counts[int(v)] += 1
    return _array_from_flat(counts, (n,), int64_dtype)

def unravel_index(indices, shape):
    if isinstance(indices, int):
        return tuple([(indices // _prod(shape[i+1:])) % shape[i] for i in range(len(shape))])
    indices = array(indices)
    result = [[] for _ in range(len(shape))]
    for flat in indices._data:
        multi = _unravel(int(flat), shape)
        for dim, idx in enumerate(multi):
            result[dim].append(idx)
    return tuple(_array_from_flat(r, (len(r),), int64_dtype) for r in result)

def ravel_multi_index(multi_index, dims, mode='raise', order='C'):
    if isinstance(multi_index, ndarray):
        multi_index = tuple(multi_index)
    result = []
    n = len(multi_index[0]._data) if isinstance(multi_index[0], ndarray) else 1
    for i in range(n):
        flat = 0
        for idx_arr, dim in zip(multi_index, dims):
            v = int(idx_arr._data[i]) if isinstance(idx_arr, ndarray) else int(idx_arr)
            flat = flat * dim + v
        result.append(flat)
    return _array_from_flat(result, (len(result),), int64_dtype)

def lexsort(keys):
    keys = [array(k)._data for k in keys]
    n = len(keys[0])
    indices = sorted(range(n), key=lambda i: tuple(k[i] for k in keys))
    return _array_from_flat(indices, (n,), int64_dtype)

def sort(a, axis=-1, kind=None, order=None):
    a = array(a).copy()
    a.sort(axis, kind, order)
    return a

def argsort(a, axis=-1, kind=None, order=None):
    return array(a).argsort(axis, kind, order)

def searchsorted(a, v, side='left', sorter=None):
    a = array(a)
    if sorter is not None:
        data = [a._data[int(i)] for i in array(sorter)._data]
    else:
        data = sorted(a._data)
    target = float(v) if not isinstance(v, ndarray) else None
    if target is not None:
        lo, hi = 0, len(data)
        while lo < hi:
            mid = (lo + hi) // 2
            if (data[mid] < target if side == 'left' else data[mid] <= target):
                lo = mid + 1
            else:
                hi = mid
        return lo
    return _array_from_flat([searchsorted(a, vi, side, sorter) for vi in v._data],
                            v._shape, int64_dtype)

def unique(ar, return_index=False, return_inverse=False, return_counts=False, axis=None):
    ar = array(ar).flatten()
    seen = []
    seen_set = set()
    indices = []
    for i, v in enumerate(ar._data):
        key = v if not isinstance(v, float) or not math.isnan(v) else 'nan'
        if key not in seen_set:
            seen_set.add(key)
            seen.append(v)
            indices.append(i)
    unique_arr = _array_from_flat(seen, (len(seen),), ar._dtype)
    result = [unique_arr]
    if return_index:
        result.append(_array_from_flat(indices, (len(indices),), int64_dtype))
    if return_inverse:
        inv = [seen.index(v) for v in ar._data]
        result.append(_array_from_flat(inv, ar._shape, int64_dtype))
    if return_counts:
        counts = [ar._data.count(v) for v in seen]
        result.append(_array_from_flat(counts, (len(counts),), int64_dtype))
    return tuple(result) if len(result) > 1 else result[0]

def pad(array_, pad_width, mode='constant', **kwargs):
    a = array(array_)
    constant_values = kwargs.get('constant_values', 0)
    if isinstance(pad_width, int):
        pad_width = [(pad_width, pad_width)] * a.ndim
    if a.ndim == 1:
        before, after = pad_width[0]
        cv = constant_values
        data = [cv] * before + list(a._data) + [cv] * after
        return _array_from_flat(data, (len(data),), a._dtype)
    raise NotImplementedError('pad for ndim > 1')

def roll(a, shift, axis=None):
    a = array(a)
    if axis is None:
        n = a.size
        shift = shift % n
        data = a._data[-shift:] + a._data[:-shift]
        return _array_from_flat(data, a._shape, a._dtype)
    raise NotImplementedError('roll with axis')

def flip(m, axis=None):
    m = array(m)
    if axis is None:
        return _array_from_flat(list(reversed(m._data)), m._shape, m._dtype)
    raise NotImplementedError

def fliplr(m):
    m = array(m)
    cols = m._shape[1]
    data = []
    for r in range(m._shape[0]):
        row = m._data[r * cols:(r + 1) * cols]
        data.extend(reversed(row))
    return _array_from_flat(data, m._shape, m._dtype)

def flipud(m):
    m = array(m)
    rows, cols = m._shape[0], _prod(m._shape[1:])
    data = []
    for r in reversed(range(rows)):
        data.extend(m._data[r * cols:(r + 1) * cols])
    return _array_from_flat(data, m._shape, m._dtype)

def rot90(m, k=1, axes=(0, 1)):
    m = array(m)
    k = k % 4
    for _ in range(k):
        m = fliplr(m).transpose()
    return m

def expand_dims(a, axis):
    a = array(a)
    sh = list(a._shape)
    sh.insert(axis, 1)
    return a.reshape(*sh)

def squeeze(a, axis=None):
    return array(a).squeeze(axis)

def broadcast_to(array_, shape):
    a = array(array_)
    out_shape = tuple(shape)
    n = _prod(out_shape)
    data = [a._data[_broadcast_index(i, a._shape, out_shape)] for i in range(n)]
    return _array_from_flat(data, out_shape, a._dtype)

def broadcast_arrays(*args, subok=False):
    arrays = [array(a) for a in args]
    shape = functools.reduce(_broadcast_shapes, [a._shape for a in arrays], ())
    return [broadcast_to(a, shape) for a in arrays]

def meshgrid(*xi, indexing='xy', sparse=False, copy=True):
    xi = [array(x).flatten() for x in xi]
    if indexing == 'xy' and len(xi) > 1:
        xi[0], xi[1] = xi[1], xi[0]
    grids = []
    for i, x in enumerate(xi):
        shape = [1] * len(xi)
        shape[i] = x.size
        grids.append(broadcast_to(x.reshape(*shape),
                                  tuple(xj.size for xj in xi)))
    return grids

def may_share_memory(a, b, max_work=None):
    return False

def shares_memory(a, b, max_work=None):
    return False

def can_cast(from_, to, casting='safe'):
    return True

def promote_types(type1, type2):
    return dtype('d')

def result_type(*arrays_and_dtypes):
    return dtype('d')

def min_scalar_type(a):
    return dtype('q')

def interp(x, xp, fp, left=None, right=None, period=None):
    x = array(x)
    xp = array(xp)._data
    fp = array(fp)._data
    def _interp_single(xi):
        if xi <= xp[0]: return left if left is not None else fp[0]
        if xi >= xp[-1]: return right if right is not None else fp[-1]
        for i in range(len(xp) - 1):
            if xp[i] <= xi <= xp[i + 1]:
                t = (xi - xp[i]) / (xp[i + 1] - xp[i])
                return fp[i] + t * (fp[i + 1] - fp[i])
        return fp[-1]
    data = [_interp_single(xi) for xi in x._data]
    return _array_from_flat(data, x._shape, float64_dtype)

def correlate(a, v, mode='valid'):
    a, v = array(a)._data, array(v)._data
    n, m = len(a), len(v)
    if mode == 'full':
        out_len = n + m - 1
        def _get(seq, i): return seq[i] if 0 <= i < len(seq) else 0
        result = [sum(_get(a, i) * _get(v, i - k + m - 1) for i in range(n + m - 1))
                  for k in range(out_len)]
    else:
        result = [sum(a[i] * v[j] for j, i in enumerate(range(k, k + m)))
                  for k in range(n - m + 1)]
    return _array_from_flat(result, (len(result),), float64_dtype)

def packbits(myarray, axis=None, bitorder='big'):
    data = array(myarray)._data
    result = []
    for i in range(0, len(data), 8):
        chunk = data[i:i+8]
        byte = 0
        for j, bit in enumerate(chunk):
            if bitorder == 'big':
                byte |= (int(bool(bit)) << (7 - j))
            else:
                byte |= (int(bool(bit)) << j)
        result.append(byte)
    return _array_from_flat(result, (len(result),), uint8_dtype)

def unpackbits(myarray, axis=None, count=None, bitorder='big'):
    data = array(myarray)._data
    result = []
    for byte in data:
        byte = int(byte)
        bits = [(byte >> (7 - i)) & 1 if bitorder == 'big' else (byte >> i) & 1
                for i in range(8)]
        result.extend(bits)
    if count is not None:
        result = result[:count]
    return _array_from_flat(result, (len(result),), uint8_dtype)

# ===========================================================================
# Section 6 – ufunc class and instances
# ===========================================================================
class ufunc:
    _nin: int
    _nout: int
    _name: str
    _fn: callable
    _identity = None

    def __init__(self, name, nin, nout, fn, identity=None, types=None):
        self._name = name
        self.__name__ = name
        self.__qualname__ = name
        self._nin = nin
        self._nout = nout
        self._fn = fn
        self._identity = identity
        self._types = types or []

    @property
    def nin(self): return self._nin
    @property
    def nout(self): return self._nout
    @property
    def nargs(self): return self._nin + self._nout
    @property
    def ntypes(self): return len(self._types) or 1
    @property
    def identity(self): return self._identity
    @property
    def types(self): return self._types
    @property
    def signature(self): return None

    def resolve_dtypes(self, *a, **kw): return None
    def _resolve_dtypes_and_context(self, *a, **kw): return None
    def _get_strided_loop(self, *a, **kw): return None

    def __repr__(self):
        return f'<ufunc {self._name!r}>'

    def __call__(self, *args, out=None, where=True, casting='same_kind',
                 order='K', dtype=None, subok=True, **kw):
        arrays = [array(a) if not isinstance(a, ndarray) else a for a in args[:self._nin]]
        if all(a._shape == () for a in arrays):
            vals = [a._data[0] for a in arrays]
            try:
                res = self._fn(*vals)
            except (ZeroDivisionError, OverflowError, ValueError):
                res = float('nan')
            if out is not None:
                out[0]._data[0] = res
                return out[0]
            if self._nout > 1:
                return tuple(_array_from_flat([r], (), float64_dtype) for r in res)
            return res
        out_shape = functools.reduce(_broadcast_shapes,
                                     [a._shape for a in arrays], ())
        n = _prod(out_shape)
        if self._nout == 1:
            result = [None] * n
            for i in range(n):
                vals = []
                for a in arrays:
                    idx = _broadcast_index(i, a._shape, out_shape) if a._shape else 0
                    vals.append(a._data[idx] if a._shape else a._data[0])
                try:
                    result[i] = self._fn(*vals)
                except (ZeroDivisionError, OverflowError):
                    result[i] = float('nan') if isinstance(vals[0], float) else 0
                except (ValueError, TypeError):
                    result[i] = float('nan')
            dt = arrays[0]._dtype if arrays else float64_dtype
            r = _array_from_flat(result, out_shape, dt)
            if out is not None:
                out[0]._data[:] = r._data
                return out[0]
            return r
        else:
            results = [[] for _ in range(self._nout)]
            for i in range(n):
                vals = []
                for a in arrays:
                    idx = _broadcast_index(i, a._shape, out_shape) if a._shape else 0
                    vals.append(a._data[idx] if a._shape else a._data[0])
                try:
                    res = self._fn(*vals)
                except Exception:
                    res = (float('nan'),) * self._nout
                if not isinstance(res, tuple):
                    res = (res,)
                for k, v in enumerate(res):
                    results[k].append(v)
            dt = arrays[0]._dtype if arrays else float64_dtype
            return tuple(_array_from_flat(r, out_shape, dt) for r in results)

    def reduce(self, a, axis=0, dtype=None, out=None, keepdims=False,
               initial=None, where=True):
        a = array(a)
        if a.ndim == 1 or axis is None:
            data = a._data
            acc = initial if initial is not None else (self._identity if self._identity is not None else data[0])
            start = 0 if initial is not None else 1
            for x in data[start:]:
                acc = self._fn(acc, x)
            return acc
        return _reduce_axis(a, lambda d: functools.reduce(self._fn, d,
                                initial if initial is not None else d[0]), axis, keepdims)

    def reduceat(self, a, indices, axis=0, dtype=None, out=None):
        return self.reduce(a, axis=axis, dtype=dtype, out=out)

    def accumulate(self, a, axis=0, dtype=None, out=None):
        a = array(a)
        if a.ndim == 1:
            acc = a._data[0]
            result = [acc]
            for x in a._data[1:]:
                acc = self._fn(acc, x)
                result.append(acc)
            return _array_from_flat(result, (len(result),), a._dtype)
        raise NotImplementedError

    def outer(self, A, B):
        A, B = array(A).flatten(), array(B).flatten()
        data = [self._fn(x, y) for x in A._data for y in B._data]
        return _array_from_flat(data, (A.size, B.size), A._dtype)

    def at(self, a, indices, b=None):
        a = array(a)
        if b is None:
            for i in array(indices)._data:
                a._data[int(i)] = self._fn(a._data[int(i)])
        else:
            b = array(b)
            for k, i in enumerate(array(indices)._data):
                bv = b._data[k] if k < b.size else b._data[0]
                a._data[int(i)] = self._fn(a._data[int(i)], bv)

def _safe_div(a, b):
    try:
        return a / b
    except ZeroDivisionError:
        return float('inf') if a >= 0 else float('-inf')

def _safe_floordiv(a, b):
    try:
        return a // b
    except ZeroDivisionError:
        return 0

def _safe_mod(a, b):
    try:
        return a % b
    except ZeroDivisionError:
        return 0

def _safe_log(x):
    try:
        return math.log(x) if x > 0 else float('-inf') if x == 0 else float('nan')
    except Exception:
        return float('nan')

def _safe_sqrt(x):
    try:
        return math.sqrt(x) if x >= 0 else float('nan')
    except Exception:
        return float('nan')

def _safe_arcsin(x):
    try:
        return math.asin(x)
    except ValueError:
        return float('nan')

def _safe_arccos(x):
    try:
        return math.acos(x)
    except ValueError:
        return float('nan')

def _safe_arctanh(x):
    try:
        return math.atanh(x)
    except ValueError:
        return float('nan')

def _safe_pow(a, b):
    try:
        return a ** b
    except (ZeroDivisionError, OverflowError):
        return float('nan')

def _logaddexp(x, y):
    if x == float('inf') and y == float('inf'):
        return float('inf')
    if x > y:
        return x + math.log1p(math.exp(y - x))
    return y + math.log1p(math.exp(x - y))

def _logaddexp2(x, y):
    return _logaddexp(x, y) / math.log(2)

def _heaviside(x, h0):
    if x < 0: return 0.0
    if x > 0: return 1.0
    return h0

def _gcd(a, b):
    a, b = int(abs(a)), int(abs(b))
    while b:
        a, b = b, a % b
    return a

def _lcm(a, b):
    g = _gcd(a, b)
    return 0 if g == 0 else abs(int(a) * int(b)) // g

def _frexp(x):
    return math.frexp(x)

def _ldexp(x, i):
    return math.ldexp(x, int(i))

def _modf(x):
    return math.modf(x)

def _divmod_ufunc(x, y):
    return divmod(x, y) if y != 0 else (float('nan'), float('nan'))

def _nextafter(x, y):
    if x == y:
        return y
    packed = _struct.pack('d', x)
    i = _struct.unpack('Q', packed)[0]
    if x < y:
        i += 1
    else:
        i -= 1
    return _struct.unpack('d', _struct.pack('Q', i))[0]

def _spacing(x):
    return _nextafter(x, float('inf')) - x

def _bitwise_count(x):
    return bin(abs(int(x))).count('1')

def _copysign(x, y):
    return math.copysign(x, y)

add = ufunc('add', 2, 1, operator.add, identity=0)
subtract = ufunc('subtract', 2, 1, operator.sub, identity=0)
multiply = ufunc('multiply', 2, 1, operator.mul, identity=1)
divide = ufunc('divide', 2, 1, _safe_div)
true_divide = ufunc('true_divide', 2, 1, _safe_div)
floor_divide = ufunc('floor_divide', 2, 1, _safe_floordiv)
power = ufunc('power', 2, 1, _safe_pow, identity=1)
float_power = ufunc('float_power', 2, 1, _safe_pow)
mod = ufunc('mod', 2, 1, _safe_mod)
remainder = ufunc('remainder', 2, 1, _safe_mod)
divmod_uf = ufunc('divmod', 2, 2, _divmod_ufunc)
negative = ufunc('negative', 1, 1, operator.neg)
positive = ufunc('positive', 1, 1, operator.pos)
absolute = ufunc('absolute', 1, 1, abs)
sign = ufunc('sign', 1, 1, lambda x: 0 if x == 0 else (1 if x > 0 else -1))
reciprocal = ufunc('reciprocal', 1, 1, lambda x: _safe_div(1.0, x))
rint = ufunc('rint', 1, 1, round)
square = ufunc('square', 1, 1, lambda x: x * x)
sqrt = ufunc('sqrt', 1, 1, _safe_sqrt)
cbrt = ufunc('cbrt', 1, 1, lambda x: math.copysign(abs(x) ** (1/3), x))
exp = ufunc('exp', 1, 1, lambda x: math.exp(x) if not isinstance(x, complex) else cmath.exp(x))
exp2 = ufunc('exp2', 1, 1, lambda x: 2.0 ** x)
expm1 = ufunc('expm1', 1, 1, math.expm1)
log = ufunc('log', 1, 1, _safe_log)
log2 = ufunc('log2', 1, 1, lambda x: _safe_log(x) / math.log(2))
log10 = ufunc('log10', 1, 1, lambda x: math.log10(x) if x > 0 else float('-inf') if x == 0 else float('nan'))
log1p = ufunc('log1p', 1, 1, lambda x: math.log1p(x) if x > -1 else float('nan'))
logaddexp = ufunc('logaddexp', 2, 1, _logaddexp)
logaddexp2 = ufunc('logaddexp2', 2, 1, _logaddexp2)
sin = ufunc('sin', 1, 1, math.sin)
cos = ufunc('cos', 1, 1, math.cos)
tan = ufunc('tan', 1, 1, math.tan)
arcsin = ufunc('arcsin', 1, 1, _safe_arcsin)
arccos = ufunc('arccos', 1, 1, _safe_arccos)
arctan = ufunc('arctan', 1, 1, math.atan)
arctan2 = ufunc('arctan2', 2, 1, math.atan2)
sinh = ufunc('sinh', 1, 1, math.sinh)
cosh = ufunc('cosh', 1, 1, math.cosh)
tanh = ufunc('tanh', 1, 1, math.tanh)
arcsinh = ufunc('arcsinh', 1, 1, math.asinh)
arccosh = ufunc('arccosh', 1, 1, lambda x: math.acosh(x) if x >= 1 else float('nan'))
arctanh = ufunc('arctanh', 1, 1, _safe_arctanh)
deg2rad = ufunc('deg2rad', 1, 1, math.radians)
rad2deg = ufunc('rad2deg', 1, 1, math.degrees)
radians = ufunc('radians', 1, 1, math.radians)
degrees = ufunc('degrees', 1, 1, math.degrees)
hypot = ufunc('hypot', 2, 1, math.hypot)
fabs = ufunc('fabs', 1, 1, lambda x: float(abs(x)))
floor = ufunc('floor', 1, 1, math.floor)
ceil = ufunc('ceil', 1, 1, math.ceil)
trunc = ufunc('trunc', 1, 1, math.trunc)
fmod = ufunc('fmod', 2, 1, math.fmod)
frexp = ufunc('frexp', 1, 2, _frexp)
ldexp = ufunc('ldexp', 2, 1, _ldexp)
modf = ufunc('modf', 1, 2, _modf)
nextafter = ufunc('nextafter', 2, 1, _nextafter)
spacing = ufunc('spacing', 1, 1, _spacing)
copysign = ufunc('copysign', 2, 1, _copysign)
heaviside = ufunc('heaviside', 2, 1, _heaviside)
isfinite = ufunc('isfinite', 1, 1, lambda x: math.isfinite(float(x)) if not isinstance(x, complex) else cmath.isfinite(x))
isinf = ufunc('isinf', 1, 1, lambda x: math.isinf(float(x)) if not isinstance(x, complex) else cmath.isinf(x))
isnan = ufunc('isnan', 1, 1, lambda x: math.isnan(float(x)) if not isinstance(x, complex) else cmath.isnan(x))
isnat = ufunc('isnat', 1, 1, lambda x: False)
signbit = ufunc('signbit', 1, 1, lambda x: math.copysign(1, x) < 0)
equal = ufunc('equal', 2, 1, operator.eq)
not_equal = ufunc('not_equal', 2, 1, operator.ne)
less = ufunc('less', 2, 1, operator.lt)
less_equal = ufunc('less_equal', 2, 1, operator.le)
greater = ufunc('greater', 2, 1, operator.gt)
greater_equal = ufunc('greater_equal', 2, 1, operator.ge)
logical_and = ufunc('logical_and', 2, 1, lambda a, b: bool(a) and bool(b), identity=True)
logical_or = ufunc('logical_or', 2, 1, lambda a, b: bool(a) or bool(b), identity=False)
logical_xor = ufunc('logical_xor', 2, 1, lambda a, b: bool(a) ^ bool(b))
logical_not = ufunc('logical_not', 1, 1, lambda a: not bool(a))
maximum = ufunc('maximum', 2, 1, max)
minimum = ufunc('minimum', 2, 1, min)
fmax = ufunc('fmax', 2, 1, lambda a, b: b if (isinstance(a, float) and math.isnan(a)) else a if (isinstance(b, float) and math.isnan(b)) else max(a, b))
fmin = ufunc('fmin', 2, 1, lambda a, b: b if (isinstance(a, float) and math.isnan(a)) else a if (isinstance(b, float) and math.isnan(b)) else min(a, b))
bitwise_and = ufunc('bitwise_and', 2, 1, operator.and_)
bitwise_or = ufunc('bitwise_or', 2, 1, operator.or_)
bitwise_xor = ufunc('bitwise_xor', 2, 1, operator.xor)
bitwise_count = ufunc('bitwise_count', 1, 1, _bitwise_count)
invert = ufunc('invert', 1, 1, operator.invert)
left_shift = ufunc('left_shift', 2, 1, operator.lshift)
right_shift = ufunc('right_shift', 2, 1, operator.rshift)
gcd = ufunc('gcd', 2, 1, _gcd, identity=0)
lcm = ufunc('lcm', 2, 1, _lcm)
conj = ufunc('conj', 1, 1, lambda x: x.conjugate() if isinstance(x, complex) else x)
conjugate = conj
clip_uf = ufunc('clip', 3, 1, lambda a, lo, hi: max(lo, min(hi, a)))
_ones_like = ufunc('_ones_like', 1, 1, lambda x: type(x)(1) if callable(type(x)) else 1)
matmul_uf = ufunc('matmul', 2, 1, lambda a, b: a * b)
vecdot = ufunc('vecdot', 2, 1, lambda a, b: a * b)
matvec = ufunc('matvec', 2, 1, lambda a, b: a * b)
vecmat = ufunc('vecmat', 2, 1, lambda a, b: a * b)

def _str_ufunc(name, nin, nout, fn):
    return ufunc(name, nin, nout, fn)

str_len = _str_ufunc('str_len', 1, 1, lambda s: len(s) if isinstance(s, (str, bytes)) else 0)
count_uf = _str_ufunc('count', 3, 1, lambda s, sub, _: s.count(sub))
find_uf = _str_ufunc('find', 4, 1, lambda s, sub, st, en: s.find(sub, int(st), int(en) if en else None))
rfind_uf = _str_ufunc('rfind', 4, 1, lambda s, sub, st, en: s.rfind(sub, int(st), int(en) if en else None))
index_uf = _str_ufunc('index', 4, 1, lambda s, sub, st, en: s.index(sub, int(st), int(en) if en else None))
rindex_uf = _str_ufunc('rindex', 4, 1, lambda s, sub, st, en: s.rindex(sub, int(st), int(en) if en else None))
startswith_uf = _str_ufunc('startswith', 4, 1, lambda s, pre, st, en: s.startswith(pre, int(st), int(en) if en else None))
endswith_uf = _str_ufunc('endswith', 4, 1, lambda s, suf, st, en: s.endswith(suf, int(st), int(en) if en else None))
isalnum_uf = _str_ufunc('isalnum', 1, 1, lambda s: s.isalnum())
isalpha_uf = _str_ufunc('isalpha', 1, 1, lambda s: s.isalpha())
isdigit_uf = _str_ufunc('isdigit', 1, 1, lambda s: s.isdigit())
isdecimal_uf = _str_ufunc('isdecimal', 1, 1, lambda s: s.isdecimal())
isnumeric_uf = _str_ufunc('isnumeric', 1, 1, lambda s: s.isnumeric() if isinstance(s, str) else False)
isspace_uf = _str_ufunc('isspace', 1, 1, lambda s: s.isspace())
islower_uf = _str_ufunc('islower', 1, 1, lambda s: s.islower())
isupper_uf = _str_ufunc('isupper', 1, 1, lambda s: s.isupper())
istitle_uf = _str_ufunc('istitle', 1, 1, lambda s: s.istitle())
_lstrip_whitespace = _str_ufunc('_lstrip_whitespace', 1, 1, lambda s: s.lstrip())
_rstrip_whitespace = _str_ufunc('_rstrip_whitespace', 1, 1, lambda s: s.rstrip())
_strip_whitespace = _str_ufunc('_strip_whitespace', 1, 1, lambda s: s.strip())
_lstrip_chars = _str_ufunc('_lstrip_chars', 2, 1, lambda s, c: s.lstrip(c))
_rstrip_chars = _str_ufunc('_rstrip_chars', 2, 1, lambda s, c: s.rstrip(c))
_strip_chars = _str_ufunc('_strip_chars', 2, 1, lambda s, c: s.strip(c))
_center = _str_ufunc('_center', 3, 1, lambda s, w, f: s.center(w, f))
_ljust = _str_ufunc('_ljust', 3, 1, lambda s, w, f: s.ljust(w, f))
_rjust = _str_ufunc('_rjust', 3, 1, lambda s, w, f: s.rjust(w, f))
_zfill = _str_ufunc('_zfill', 2, 1, lambda s, w: s.zfill(w))
_expandtabs = _str_ufunc('_expandtabs', 2, 1, lambda s, t: s.expandtabs(t))
_expandtabs_length = _str_ufunc('_expandtabs_length', 2, 1, lambda s, t: len(s.expandtabs(t)))
_slice = _str_ufunc('_slice', 4, 1, lambda s, st, en, step: s[int(st):int(en):int(step)])
_replace = _str_ufunc('_replace', 4, 1, lambda s, old, new, n: s.replace(old, new, int(n)))
_partition = _str_ufunc('_partition', 2, 1, lambda s, sep: s.partition(sep))
_rpartition = _str_ufunc('_rpartition', 2, 1, lambda s, sep: s.rpartition(sep))
_partition_index = _str_ufunc('_partition_index', 3, 1, lambda s, sep, _: s.partition(sep))
_rpartition_index = _str_ufunc('_rpartition_index', 3, 1, lambda s, sep, _: s.rpartition(sep))
_arg = _str_ufunc('_arg', 1, 1, lambda s: s)

# ===========================================================================
# Section 7 – Extra math functions
# ===========================================================================
def sum_(a, axis=None, dtype=None, out=None, keepdims=False, initial=None, where=True):
    a = array(a)
    if axis is None:
        s = sum(a._data)
        return s if initial is None else s + initial
    return _reduce_axis(a, sum, axis, keepdims)

def cumsum(a, axis=None, dtype=None, out=None):
    a = array(a)
    data = a._data
    acc = 0
    result = []
    for x in data:
        acc += x
        result.append(acc)
    return _array_from_flat(result, a._shape if axis is not None else (len(result),), float64_dtype)

def cumprod(a, axis=None, dtype=None, out=None):
    a = array(a)
    acc = 1
    result = []
    for x in a._data:
        acc *= x
        result.append(acc)
    return _array_from_flat(result, (len(result),), float64_dtype)

def diff(a, n=1, axis=-1, prepend=None, append=None):
    a = array(a)
    for _ in range(n):
        data = [a._data[i+1] - a._data[i] for i in range(len(a._data) - 1)]
        a = _array_from_flat(data, (len(data),), float64_dtype)
    return a

def gradient(f, *varargs, **kwargs):
    f = array(f)
    if f.ndim == 1:
        n = f.size
        dx = varargs[0] if varargs else 1.0
        g = [0.0] * n
        if n > 1:
            g[0] = (f._data[1] - f._data[0]) / dx
            g[-1] = (f._data[-1] - f._data[-2]) / dx
            for i in range(1, n - 1):
                g[i] = (f._data[i + 1] - f._data[i - 1]) / (2 * dx)
        return _array_from_flat(g, (n,), float64_dtype)
    raise NotImplementedError

def clip_fn(a, a_min=None, a_max=None, out=None, **kwargs):
    return array(a).clip(a_min, a_max)

def around(a, decimals=0, out=None):
    a = array(a)
    data = [round(x, decimals) for x in a._data]
    return _array_from_flat(data, a._shape, a._dtype)

def round_(a, decimals=0, out=None):
    return around(a, decimals)

def fix(x, out=None):
    x = array(x)
    data = [math.trunc(v) for v in x._data]
    return _array_from_flat(data, x._shape, float64_dtype)

def angle(z, deg=False):
    z = array(z)
    data = [math.degrees(cmath.phase(v)) if deg else cmath.phase(v) for v in z._data]
    return _array_from_flat(data, z._shape, float64_dtype)

def real(a):
    return array(a).real

def imag(a):
    return array(a).imag

def abs_(a):
    a = array(a)
    return _array_from_flat([abs(x) for x in a._data], a._shape, a._dtype)

def choose(a, choices, out=None, mode='raise'):
    a = array(a)
    choices = [array(c) for c in choices]
    data = [choices[int(i)]._data[0] if choices[int(i)].size == 1 else
            choices[int(i)]._data[j] for j, i in enumerate(a._data)]
    return _array_from_flat(data, a._shape, choices[0]._dtype)

def piecewise(x, condlist, funclist, *args, **kw):
    x = array(x)
    result = [0.0] * x.size
    for cond, fn in zip(condlist, funclist):
        cond = array(cond)
        for i, c in enumerate(cond._data):
            if c:
                result[i] = fn(x._data[i]) if callable(fn) else fn
    return _array_from_flat(result, x._shape, float64_dtype)

def select(condlist, choicelist, default=0):
    condlist = [array(c) for c in condlist]
    choicelist = [array(ch) for ch in choicelist]
    n = condlist[0].size
    result = [default] * n
    for i in range(n):
        for cond, choice in zip(condlist, choicelist):
            if cond._data[i]:
                result[i] = choice._data[i]
                break
    return _array_from_flat(result, condlist[0]._shape, float64_dtype)

def apply_along_axis(func1d, axis, arr, *args, **kwargs):
    arr = array(arr)
    sub = _prod(arr._shape[axis + 1:]) if axis + 1 < arr.ndim else 1
    outer = _prod(arr._shape[:axis])
    n_axis = arr._shape[axis]
    result = []
    for o in range(outer):
        for s in range(sub):
            chunk = [arr._data[o * n_axis * sub + j * sub + s] for j in range(n_axis)]
            sub_arr = _array_from_flat(chunk, (n_axis,), arr._dtype)
            out = func1d(sub_arr, *args, **kwargs)
            if isinstance(out, ndarray):
                result.extend(out._data)
            else:
                result.append(out)
    return _array_from_flat(result, (len(result),), float64_dtype)

def place(arr, mask, vals):
    arr = array(arr)
    vals = list(array(vals)._data)
    vi = 0
    for i, m in enumerate(array(mask)._data):
        if m:
            arr._data[i] = vals[vi % len(vals)]
            vi += 1

def putmask(a, mask, values):
    a = array(a)
    mask = array(mask)
    values = array(values)
    vi = 0
    for i, m in enumerate(mask._data):
        if m:
            a._data[i] = values._data[vi % values.size]
            vi += 1

def extract(condition, arr):
    arr = array(arr)
    condition = array(condition)
    data = [v for c, v in zip(condition._data, arr._data) if c]
    return _array_from_flat(data, (len(data),), arr._dtype)

def in1d(ar1, ar2, assume_unique=False, invert=False):
    ar1, ar2 = array(ar1).flatten(), array(ar2).flatten()
    s = set(ar2._data)
    data = [(v not in s) if invert else (v in s) for v in ar1._data]
    return _array_from_flat(data, ar1._shape, bool_dtype)

def isin(element, test_elements, assume_unique=False, invert=False):
    return in1d(element, test_elements, assume_unique, invert).reshape(array(element)._shape)

def intersect1d(ar1, ar2, assume_unique=False, return_indices=False):
    a, b = set(array(ar1)._data), set(array(ar2)._data)
    result = sorted(a & b)
    return _array_from_flat(result, (len(result),), float64_dtype)

def union1d(ar1, ar2):
    a = sorted(set(array(ar1)._data) | set(array(ar2)._data))
    return _array_from_flat(a, (len(a),), float64_dtype)

def setdiff1d(ar1, ar2, assume_unique=False):
    a = sorted(set(array(ar1)._data) - set(array(ar2)._data))
    return _array_from_flat(a, (len(a),), float64_dtype)

def setxor1d(ar1, ar2, assume_unique=False):
    a = sorted(set(array(ar1)._data) ^ set(array(ar2)._data))
    return _array_from_flat(a, (len(a),), float64_dtype)

# ===========================================================================
# Section 8 – broadcast / nditer / flatiter stubs
# ===========================================================================
class broadcast:
    def __init__(self, *args):
        self._args = [array(a) for a in args]
        self._shape = functools.reduce(_broadcast_shapes,
                                       [a._shape for a in self._args], ())
        self._size = _prod(self._shape)
        self._index = 0

    @property
    def shape(self): return self._shape
    @property
    def size(self): return self._size
    @property
    def ndim(self): return len(self._shape)
    @property
    def nd(self): return len(self._shape)
    @property
    def index(self): return self._index
    @property
    def numiter(self): return len(self._args)
    @property
    def iters(self): return tuple(iter(a._data) for a in self._args)

    def reset(self):
        self._index = 0

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index >= self._size:
            raise StopIteration
        i = self._index
        self._index += 1
        result = []
        for a in self._args:
            idx = _broadcast_index(i, a._shape, self._shape) if a._shape else 0
            result.append(a._data[idx] if a._shape else a._data[0])
        return tuple(result)

class flatiter:
    def __init__(self, arr):
        self._arr = arr
        self._idx = 0

    @property
    def index(self): return self._idx
    @property
    def base(self): return self._arr
    @property
    def coords(self): return _unravel(self._idx, self._arr._shape)

    def __iter__(self): return self
    def __next__(self):
        if self._idx >= self._arr.size:
            raise StopIteration
        v = self._arr._data[self._idx]
        self._idx += 1
        return v
    def __getitem__(self, i): return self._arr._data[i]
    def __setitem__(self, i, v): self._arr._data[i] = v
    def __len__(self): return self._arr.size
    def copy(self): return list(self._arr._data)
    def __array__(self, dtype=None, copy=None): return array(list(self._arr._data), dtype=dtype)

class nditer:
    def __init__(self, op, flags=None, op_flags=None, op_dtypes=None,
                 order='K', casting='safe', op_axes=None, itershape=None,
                 buffersize=0):
        if isinstance(op, ndarray):
            op = [op]
        self._arrays = [array(o) for o in op]
        self._flags = flags or []
        self._idx = 0
        self._shape = functools.reduce(_broadcast_shapes,
                                       [a._shape for a in self._arrays], ())
        self._size = _prod(self._shape)
        self.finished = self._size == 0
        self.iterindex = 0
        self.itersize = self._size

    @property
    def shape(self): return self._shape
    @property
    def ndim(self): return len(self._shape)
    @property
    def iterrange(self): return range(self._size)

    def __iter__(self):
        self._idx = 0
        self.finished = False
        return self

    def __next__(self):
        if self._idx >= self._size:
            self.finished = True
            raise StopIteration
        i = self._idx
        self._idx += 1
        self.iterindex = i
        result = []
        for a in self._arrays:
            idx = _broadcast_index(i, a._shape, self._shape) if a._shape else 0
            v = a._data[idx] if a._shape else a._data[0]
            result.append(_array_from_flat([v], (), a._dtype))
        return result if len(result) > 1 else result[0]

    def __enter__(self): return self
    def __exit__(self, *a): pass

    def iternext(self):
        try:
            next(self)
            return True
        except StopIteration:
            return False

    @property
    def value(self):
        i = self._idx - 1
        return [_array_from_flat([a._data[_broadcast_index(i, a._shape, self._shape)]],
                                  (), a._dtype) for a in self._arrays]

    @property
    def operands(self):
        return self._arrays

    def copy(self): return self
    def remove_axis(self, i): pass
    def remove_multi_index(self): pass
    def enable_flat_iteration(self): pass
    def enable_external_loop(self): pass
    def debug_print(self): pass
    def reset(self): self._idx = 0; self.finished = False
    def close(self): pass

    def __getitem__(self, idx): return self.operands[idx]
    def __setitem__(self, idx, val): self.operands[idx]._data = array(val)._data

class nested_iters:
    def __init__(self, op, axes, flags=None, op_flags=None,
                 op_dtypes=None, order='K', casting='safe', buffersize=0):
        self._iters = [nditer(op, flags, op_flags) for _ in axes]
    def __iter__(self): return iter(self._iters)
    def __getitem__(self, i): return self._iters[i]
    def close(self): pass

class _array_converter:
    def __init__(self, *args):
        self._args = [array(a) for a in args]
        self._types = [type(a) for a in args]
    def __call__(self, arr, **kw):
        return array(arr)
    @property
    def scalar_input(self):
        return all(a.ndim == 0 for a in self._args)
    def as_arrays(self, *args, **kw):
        return tuple(array(a) for a in (args if args else self._args))
    def result_type(self, *args, ensure_inexact=False, **kw):
        dt = self._args[0]._dtype if self._args else dtype('float64')
        if ensure_inexact and dt.kind in ('i', 'u', 'b'):
            dt = dtype('d')
        return dt
    def wrap(self, arr, like=None, **kw):
        return array(arr)

class _ArrayFunctionDispatcher:
    def __init__(self, dispatcher, func=None):
        self._dispatcher = dispatcher
        self._func = func
    def __call__(self, *args, **kw):
        if self._func is not None:
            return self._func(*args, **kw)
        return self._dispatcher(*args, **kw)

# ===========================================================================
# Section 9 – datetime helpers (stubs)
# ===========================================================================
class busdaycalendar:
    weekmask = 'Mon Tue Wed Thu Fri'
    holidays = None
    def __init__(self, weekmask='Mon Tue Wed Thu Fri', holidays=None):
        self.weekmask = weekmask
        self.holidays = array(holidays or [])

def busday_offset(dates, offsets, roll='raise', weekmask='Mon Tue Wed Thu Fri',
                  holidays=None, busdaycal=None, out=None):
    return dates

def busday_count(begindates, enddates, weekmask='Mon Tue Wed Thu Fri',
                 holidays=None, busdaycal=None, out=None):
    return 0

def is_busday(dates, weekmask='Mon Tue Wed Thu Fri', holidays=None,
              busdaycal=None, out=None):
    return True

def datetime_data(dtype):
    return ('D', 1)

def datetime_as_string(arr, unit=None, timezone='naive', casting='same_kind'):
    return array([str(v) for v in array(arr)._data])

# ===========================================================================
# Section 10 – typeinfo dict
# ===========================================================================
class _TypeInfo:
    __slots__ = ('type', 'num', 'bits', 'char', 'str', 'kind',
                 'flags', 'alignment', 'max', 'min')
    def __init__(self, tp, num, bits, char, str_, kind,
                 flags=0, alignment=1, max_val=None, min_val=None):
        self.type = tp
        self.num = num
        self.bits = bits
        self.char = char
        self.str = str_
        self.kind = kind
        self.flags = flags
        self.alignment = alignment
        self.max = max_val
        self.min = min_val

    def __len__(self):
        return 6 if self.max is not None else 4

    def __repr__(self):
        return f'_TypeInfo({self.type.__name__}, bits={self.bits})'

def _make_typeinfo():
    mapping = [
        ('bool', bool_, 0, '?', 'b', 'NPY_BOOL'),
        ('int8', int8, 1, 'b', 'i', 'NPY_BYTE'),
        ('uint8', uint8, 2, 'B', 'u', 'NPY_UBYTE'),
        ('int16', int16, 3, 'h', 'i', 'NPY_SHORT'),
        ('uint16', uint16, 4, 'H', 'u', 'NPY_USHORT'),
        ('int32', int32, 5, 'i', 'i', 'NPY_INT'),
        ('uint32', uint32, 6, 'I', 'u', 'NPY_UINT'),
        ('int64', int64, 7, 'q', 'i', 'NPY_LONG'),
        ('uint64', uint64, 8, 'Q', 'u', 'NPY_ULONG'),
        ('float16', float16, 23, 'e', 'f', 'NPY_HALF'),
        ('float32', float32, 11, 'f', 'f', 'NPY_FLOAT'),
        ('float64', float64, 12, 'd', 'f', 'NPY_DOUBLE'),
        ('longdouble', float64, 13, 'g', 'f', 'NPY_LONGDOUBLE'),
        ('complex64', complex64, 14, 'F', 'c', 'NPY_CFLOAT'),
        ('complex128', complex128, 15, 'D', 'c', 'NPY_CDOUBLE'),
        ('clongdouble', complex128, 16, 'G', 'c', 'NPY_CLONGDOUBLE'),
        ('object_', object_, 17, 'O', 'O', 'NPY_OBJECT'),
        ('bytes_', bytes_, 18, 'S', 'S', 'NPY_STRING'),
        ('str_', str_, 19, 'U', 'U', 'NPY_UNICODE'),
        ('void', void, 20, 'V', 'V', 'NPY_VOID'),
        ('datetime64', datetime64, 21, 'M', 'M', 'NPY_DATETIME'),
        ('timedelta64', timedelta64, 22, 'm', 'm', 'NPY_TIMEDELTA'),
        ('longlong', int64, 7, 'q', 'i', 'NPY_LONGLONG'),
        ('ulonglong', uint64, 8, 'Q', 'u', 'NPY_ULONGLONG'),
        ('byte', int8, 1, 'b', 'i', None),
        ('ubyte', uint8, 2, 'B', 'u', None),
        ('short', int16, 3, 'h', 'i', None),
        ('ushort', uint16, 4, 'H', 'u', None),
        ('intc', int32, 5, 'i', 'i', None),
        ('uintc', uint32, 6, 'I', 'u', None),
        ('int', int32, 5, 'i', 'i', None),
        ('uint', uint32, 6, 'I', 'u', None),
        ('float', float32, 11, 'f', 'f', None),
        ('double', float64, 12, 'd', 'f', None),
        ('cfloat', complex64, 14, 'F', 'c', None),
        ('cdouble', complex128, 15, 'D', 'c', None),
        ('string', bytes_, 18, 'S', 'S', None),
        ('unicode', str_, 19, 'U', 'U', None),
        ('object', object_, 17, 'O', 'O', None),
        ('long', int64, 7, 'q', 'i', None),
        ('ulong', uint64, 8, 'Q', 'u', None),
        ('intp', int64, 7, 'q', 'i', None),
        ('uintp', uint64, 8, 'Q', 'u', None),
        ('float_', float64, 12, 'd', 'f', None),
        ('complex_', complex128, 15, 'D', 'c', None),
        ('int_', int64, 7, 'q', 'i', None),
    ]
    ti = {}
    for cls in (generic, number, integer, signedinteger, unsignedinteger,
                inexact, floating, complexfloating, flexible, character):
        ti[cls.__name__] = cls
    for name, cls, num, char, kind, npy_alias in mapping:
        dt = dtype(char)
        bits = dt.itemsize * 8
        max_val = min_val = None
        if kind == 'i':
            max_val = (1 << (bits - 1)) - 1
            min_val = -(1 << (bits - 1))
        elif kind == 'u':
            max_val = (1 << bits) - 1
            min_val = 0
        elif kind == 'b':
            max_val = 1
            min_val = 0
        alignment = dt.itemsize
        info = _TypeInfo(cls, num, bits, char, dt.str, kind,
                         alignment=alignment, max_val=max_val, min_val=min_val)
        ti[name] = info
        if npy_alias:
            ti[npy_alias] = info
    return ti

typeinfo = _make_typeinfo()

# ===========================================================================
# Section 11 – scalar() and _reconstruct()
# ===========================================================================
def scalar(dtype, obj):
    dt = _coerce_dtype(dtype)
    if isinstance(obj, bytes):
        try:
            v = _struct.unpack('=' + dt.str[1], obj[:dt.itemsize])[0]
        except Exception:
            v = 0
        return dt.type(v)
    return dt.type(obj)

def _reconstruct(subtype, basetype, state):
    obj = ndarray.__new__(ndarray)
    obj._data = []
    obj._shape = (0,)
    obj._dtype = float64_dtype
    obj._order = 'C'
    if isinstance(state, tuple) and len(state) >= 5:
        version, shape, dt, forder, data = state[:5]
        obj._shape = tuple(shape)
        obj._dtype = _coerce_dtype(dt) or float64_dtype
        obj._order = 'F' if forder else 'C'
        if isinstance(data, bytes):
            fmt = obj._dtype.str[1]
            sz = obj._dtype.itemsize
            n = len(data) // sz
            obj._data = []
            for i in range(n):
                chunk = data[i * sz:(i + 1) * sz]
                try:
                    v = _struct.unpack('=' + fmt, chunk)[0]
                except Exception:
                    v = 0
                obj._data.append(obj._dtype.type(v))
        else:
            obj._data = list(data) if hasattr(data, '__iter__') else [data]
    return obj

def _get_ndarray_c_version():
    return 0x01000009

def normalize_axis_index(axis, ndim, msg_prefix=None):
    if axis < 0:
        axis += ndim
    if not 0 <= axis < ndim:
        raise AxisError(f'axis {axis} is out of bounds for array of dimension {ndim}')
    return axis

class AxisError(ValueError, IndexError):
    pass

def dragon4_positional(x, *a, **k): return str(x)
def dragon4_scientific(x, *a, **k): return f'{x:e}'
def format_longfloat(x, precision=None): return str(x)
def _add_newdoc_ufunc(ufunc_obj, doc): pass
def add_docstring(obj, docstring): pass
def _get_extobj_dict(): return {}
def _make_extobj(**kw): return None

import contextvars as _cv
_extobj_contextvar = _cv.ContextVar("_extobj_contextvar", default=None)

def _get_implementing_args(relevant_args, array_like_args):
    return []

def _get_castingimpl(from_, to):
    return None

def _get_sfloat_dtype(): return float32_dtype

def _populate_finfo_constants(finfo_obj, dt):
    _data = {
        'e': (10, 5), 'f': (23, 8), 'd': (52, 11), 'g': (52, 11),
    }
    char = getattr(dt, 'char', getattr(dt, '_char', 'd'))
    m, e = _data.get(char, (52, 11))
    r = 2
    precision = int(m * math.log10(r)) + 1
    eps = r ** -m
    maxexp = 2 ** (e - 1)
    minexp = 2 - maxexp
    max_val = (2.0 - 2.0**(-m)) * 2.0**(maxexp - 1)
    tiny = 2.0 ** (minexp - 1)
    sub = tiny * 2.0 ** (-(m - 1))
    finfo_obj.eps = eps
    finfo_obj.epsneg = eps / r
    finfo_obj.tiny = tiny
    finfo_obj.smallest_normal = tiny
    finfo_obj.smallest_subnormal = sub
    finfo_obj.max = max_val
    finfo_obj.min = -max_val
    finfo_obj.precision = precision
    finfo_obj._radix = r
    finfo_obj.maxexp = maxexp
    finfo_obj.minexp = minexp
    finfo_obj.nexp = e
    finfo_obj.nmant = m

def _blas_supports_fpe(arg=None): return False
def _set_madvise_hugepage(enabled): pass
def _get_madvise_hugepage(): return False
def _set_numpy_warn_if_no_mem_policy(val): pass
def _reload_guard(): pass
def set_datetimeparse_function(f): pass
def set_typeDict(d): pass
def _monotonicity(arr): return 0

def compare_chararrays(a1, a2, cmp, rstrip):
    return array(a1) == array(a2)

def _vec_string(a, dt, op, args=()):
    a = array(a)
    fn = getattr(str, op, None)
    if fn is None:
        return a.copy()
    data = [fn(str(v), *args) for v in a._data]
    return _array_from_flat(data, a._shape, dtype('U'))

def _discover_array_parameters(a, **kw):
    a = array(a)
    return {'dtype': a._dtype, 'shape': a._shape}

def get_handler_name(): return 'default'
def get_handler_version(): return 1
def _unique_hash(a): return id(a)
def _load_from_filelike(fp, dtype, sep, conv): raise NotImplementedError

def fromfile(file, dtype=float, count=-1, sep='', offset=0, *, like=None):
    raise NotImplementedError('fromfile() requires a filesystem')

def from_dlpack(obj, /, *, device=None, copy=None):
    raise NotImplementedError('from_dlpack() requires DLPack support')

def correlate2(a, v, mode='full'):
    return correlate(a, v, mode)

def interp_complex(x, xp, fp, left=None, right=None):
    return interp(x, xp, fp, left, right)

def _place(a, mask, vals):
    place(a, mask, vals)

tracemalloc_domain = 390

DATETIMEUNITS = {
    'Y': 0, 'M': 1, 'W': 2, 'D': 3, 'h': 4, 'm': 5, 's': 6,
    'ms': 7, 'us': 8, 'ns': 9, 'ps': 10, 'fs': 11, 'as': 12,
}

ALLOW_THREADS = 1
BUFSIZE = 8192
CLIP = 0
WRAP = 1
RAISE = 2
MAXDIMS = 64
MAY_SHARE_BOUNDS = 0
MAY_SHARE_EXACT = -1
FLOATING_POINT_SUPPORT = 1
FPE_DIVIDEBYZERO = 1
FPE_INVALID = 8
FPE_OVERFLOW = 4
FPE_UNDERFLOW = 2
ITEM_HASOBJECT = 1
ITEM_IS_POINTER = 4
LIST_PICKLE = 2
NEEDS_INIT = 8
NEEDS_PYAPI = 16
USE_GETITEM = 32
USE_SETITEM = 64
UFUNC_BUFSIZE_DEFAULT = 8192
UFUNC_PYVALS_NAME = 'UFUNC_PYVALS'
NAN = float('nan')
NINF = float('-inf')
PINF = float('inf')
NZERO = -0.0
PZERO = 0.0
e = math.e
pi = math.pi
euler_gamma = 0.5772156649015328865
_ARRAY_API = 16
_UFUNC_API = 4
_flagdict = {
    'C_CONTIGUOUS': 0x1, 'F_CONTIGUOUS': 0x2,
    'OWNDATA': 0x4, 'WRITEABLE': 0x400,
    'ALIGNED': 0x200, 'WRITEBACKIFCOPY': 0x2000,
}

# ===========================================================================
# Section 11.5 – iinfo and finfo
# ===========================================================================
class finfo:
    """Floating point limits."""
    _finfo_cache = {}

    def __new__(cls, dtype):
        dt = _coerce_dtype(dtype)
        key = dt.char
        if key in cls._finfo_cache:
            return cls._finfo_cache[key]
        obj = object.__new__(cls)
        obj.dtype = dt
        obj.kind = dt.kind
        _data = {
            'e': (10, 5), 'f': (23, 8), 'd': (52, 11), 'g': (52, 11),
        }
        m, e_bits = _data.get(key, (52, 11))
        r = 2
        obj.precision = int(m * math.log10(r)) + 1
        obj.eps = r ** -m
        obj.epsneg = obj.eps / r
        obj.maxexp = 2 ** (e_bits - 1)
        obj.minexp = 2 - obj.maxexp
        obj.max = (2.0 - 2.0**(-m)) * 2.0**(obj.maxexp - 1)
        obj.min = -obj.max
        obj.tiny = 2.0 ** (obj.minexp - 1)
        obj.smallest_normal = obj.tiny
        obj.smallest_subnormal = obj.tiny * 2.0 ** (-(m - 1))
        obj._radix = r
        obj.nexp = e_bits
        obj.nmant = m
        obj.bits = dt.itemsize * 8
        cls._finfo_cache[key] = obj
        return obj

    def __repr__(self):
        return f"finfo(dtype={self.dtype.name})"


class iinfo:
    """Integer limits."""
    _iinfo_cache = {}

    # Concrete struct-format chars → abstract kind
    _SIGNED   = frozenset({'b', 'h', 'i', 'l', 'q'})   # int8..int64
    _UNSIGNED = frozenset({'B', 'H', 'I', 'L', 'Q'})   # uint8..uint64
    _BOOL     = frozenset({'?'})

    def __new__(cls, dtype):
        dt = _coerce_dtype(dtype)
        key = dt.char
        if key in cls._iinfo_cache:
            return cls._iinfo_cache[key]

        # Map concrete struct char → abstract kind
        if dt.char in cls._SIGNED:
            abstract_kind = 'i'
        elif dt.char in cls._UNSIGNED:
            abstract_kind = 'u'
        elif dt.char in cls._BOOL or dt.kind == 'b':
            abstract_kind = 'b'
        elif dt.kind in ('i', 'u'):
            abstract_kind = dt.kind
        else:
            raise ValueError(f"Invalid integer data type {dt.kind!r}.")

        obj = object.__new__(cls)
        obj.dtype = dt
        obj.kind = abstract_kind
        obj.bits = dt.itemsize * 8

        if abstract_kind == 'u':
            obj.min = 0
            obj.max = (1 << obj.bits) - 1
        elif abstract_kind == 'i':
            obj.min = -(1 << (obj.bits - 1))
            obj.max = (1 << (obj.bits - 1)) - 1
        elif abstract_kind == 'b':
            obj.min = 0
            obj.max = 1

        cls._iinfo_cache[key] = obj
        return obj

    def __repr__(self):
        return f"iinfo(dtype={self.dtype.name})"

# ===========================================================================
# Section 12 – numpy._core._multiarray_umath module object
# ===========================================================================
_MULTIARRAY_UMATH_ATTRS = {
    'ndarray': ndarray, 'dtype': dtype, 'ufunc': ufunc,
    'broadcast': broadcast, 'flatiter': flatiter, 'nditer': nditer,
    'flagsobj': flagsobj, 'generic': generic, 'number': number,
    'integer': integer, 'signedinteger': signedinteger,
    'unsignedinteger': unsignedinteger, 'inexact': inexact,
    'floating': floating, 'complexfloating': complexfloating,
    'flexible': flexible, 'character': character,
    'busdaycalendar': busdaycalendar, '_array_converter': _array_converter,
    '_ArrayFunctionDispatcher': _ArrayFunctionDispatcher,
    'StringDType': str_, 'error': Exception,
    # CPU feature detection stubs
    'cpu_features': {},
    'cpu_baseline': [],
    'cpu_dispatch': [],
    'cpu_targets_info': {},
    # Double-underscore versions (required by numpy 2.x)
    '__cpu_features__': {},
    '__cpu_baseline__': [],
    '__cpu_dispatch__': [],
    'version': '2.4.2', 'git_version': '0' * 40,
    'geterrobj': lambda: [0, 0, None], 'seterrobj': lambda obj: None,
    'DATETIMEUNITS': {u: i for i, u in enumerate(
        ['Y','M','W','D','h','m','s','ms','us','ns','ps','fs','as'])},
    'ITEM_HASOBJECT': 1, 'ITEM_IS_POINTER': 4, 'LIST_PICKLE': 2,
    'NEEDS_INIT': 8, 'NEEDS_PYAPI': 16, 'USE_GETITEM': 32, 'USE_SETITEM': 64,
    'correlate2': lambda a, v, mode='full': array([]),
    'count_nonzero': lambda a, axis=None, keepdims=False: sum(1 for x in array(a)._data if x),
    'dragon4_positional': lambda *a, **kw: '',
    'dragon4_scientific': lambda *a, **kw: '',
    'format_longfloat': lambda x: str(x),
    'interp': lambda x, xp, fp, left=None, right=None, period=None: 0.0,
    'interp_complex': lambda x, xp, fp, left=None, right=None, period=None: 0j,
    'packbits': lambda a, axis=None, bitorder='big': array([], dtype=dtype('B')),
    'unpackbits': lambda a, axis=None, count=None, bitorder='big': array([], dtype=dtype('B')),
    'ravel_multi_index': lambda multi_index, dims, mode='raise', order='C': array([]),
    'unravel_index': lambda indices, shape, order='C': tuple(array([]) for _ in shape),
    'bincount': lambda x, weights=None, minlength=0: array([]),
    'set_datetimeparse_function': lambda f=None: None,
    'set_legacy_print_mode': lambda val: None,
    'set_numeric_ops': lambda **ops: {},
    'set_string_function': lambda f, repr=True: None,
    'fastCopyAndTranspose': lambda a: array(a).T.copy(),
    '_flagdict': {'C_CONTIGUOUS':1,'F_CONTIGUOUS':2,'OWNDATA':4,
                  'WRITEABLE':8,'ALIGNED':16,'WRITEBACKIFCOPY':32},
    '_place': lambda a, m, v: None,
    '_vec_string': lambda a, dt, op, args=(): array([], dtype=dt),
    '_ARRAY_API': None, '_monotonicity': lambda a: True,
    '_get_promotion_state': lambda: 'weak', '_set_promotion_state': lambda s: None,
    '_using_numpy2_behavior': lambda: False,
    '_get_implementing_args': lambda relevant_args, array_types, args, kwargs: [],
    'tracemalloc_domain': 1, '_load_from_filelike': lambda *a, **kw: None,
    '_get_castingimpl': lambda *a, **kw: None,
    '_discover_array_parameters': lambda *a, **kw: (None, None),
    '_UFUNC_API': None, '_add_newdoc_ufunc': lambda *a, **kw: None,
    '_ones_like': ones,
    'add': add, 'subtract': subtract, 'multiply': multiply,
    'divide': divide, 'true_divide': true_divide, 'floor_divide': floor_divide,
    'power': power, 'float_power': float_power, 'mod': mod,
    'remainder': remainder, 'divmod': divmod_uf,
    'negative': negative, 'positive': positive, 'absolute': absolute,
    'sign': sign, 'reciprocal': reciprocal, 'rint': rint,
    'square': square, 'sqrt': sqrt, 'cbrt': cbrt,
    'exp': exp, 'exp2': exp2, 'expm1': expm1,
    'log': log, 'log2': log2, 'log10': log10, 'log1p': log1p,
    'logaddexp': logaddexp, 'logaddexp2': logaddexp2,
    'sin': sin, 'cos': cos, 'tan': tan,
    'arcsin': arcsin, 'arccos': arccos, 'arctan': arctan, 'arctan2': arctan2,
    'sinh': sinh, 'cosh': cosh, 'tanh': tanh,
    'arcsinh': arcsinh, 'arccosh': arccosh, 'arctanh': arctanh,
    'deg2rad': deg2rad, 'rad2deg': rad2deg, 'radians': radians, 'degrees': degrees,
    'hypot': hypot, 'fabs': fabs, 'floor': floor, 'ceil': ceil, 'trunc': trunc,
    'fmod': fmod, 'frexp': frexp, 'ldexp': ldexp, 'modf': modf,
    'nextafter': nextafter, 'spacing': spacing, 'copysign': copysign,
    'heaviside': heaviside,
    'isfinite': isfinite, 'isinf': isinf, 'isnan': isnan, 'isnat': isnat,
    'signbit': signbit,
    'equal': equal, 'not_equal': not_equal,
    'less': less, 'less_equal': less_equal,
    'greater': greater, 'greater_equal': greater_equal,
    'logical_and': logical_and, 'logical_or': logical_or,
    'logical_xor': logical_xor, 'logical_not': logical_not,
    'maximum': maximum, 'minimum': minimum, 'fmax': fmax, 'fmin': fmin,
    'bitwise_and': bitwise_and, 'bitwise_or': bitwise_or,
    'bitwise_xor': bitwise_xor, 'bitwise_count': bitwise_count,
    'invert': invert, 'left_shift': left_shift, 'right_shift': right_shift,
    'gcd': gcd, 'lcm': lcm,
    'conj': conj, 'conjugate': conjugate,
    'clip': clip_uf, '_ones_like': _ones_like,
    'matmul': matmul_uf, 'vecdot': vecdot, 'matvec': matvec, 'vecmat': vecmat,
    'str_len': str_len, 'count': count_uf, 'find': find_uf, 'rfind': rfind_uf,
    'index': index_uf, 'rindex': rindex_uf,
    'startswith': startswith_uf, 'endswith': endswith_uf,
    'isalnum': isalnum_uf, 'isalpha': isalpha_uf, 'isdigit': isdigit_uf,
    'isdecimal': isdecimal_uf, 'isnumeric': isnumeric_uf,
    'isspace': isspace_uf, 'islower': islower_uf, 'isupper': isupper_uf,
    'istitle': istitle_uf,
    '_lstrip_whitespace': _lstrip_whitespace, '_rstrip_whitespace': _rstrip_whitespace,
    '_strip_whitespace': _strip_whitespace,
    '_lstrip_chars': _lstrip_chars, '_rstrip_chars': _rstrip_chars,
    '_strip_chars': _strip_chars,
    '_center': _center, '_ljust': _ljust, '_rjust': _rjust,
    '_zfill': _zfill, '_expandtabs': _expandtabs,
    '_expandtabs_length': _expandtabs_length,
    '_slice': _slice, '_replace': _replace,
    '_partition': _partition, '_rpartition': _rpartition,
    '_partition_index': _partition_index, '_rpartition_index': _rpartition_index,
    '_arg': _arg,
    'array': array, 'zeros': zeros, 'ones': ones, 'empty': empty,
    'arange': arange, 'frombuffer': frombuffer, 'fromiter': fromiter,
    'fromstring': fromstring, 'frompyfunc': frompyfunc,
    'concatenate': concatenate, 'dot': dot, 'inner': inner,
    'vdot': vdot, 'matmul': matmul, 'c_einsum': c_einsum,
    'count_nonzero': count_nonzero, 'bincount': bincount,
    'can_cast': can_cast, 'promote_types': promote_types,
    'result_type': result_type, 'min_scalar_type': min_scalar_type,
    'interp': interp, 'interp_complex': interp_complex,
    'correlate': correlate, 'correlate2': correlate2,
    'lexsort': lexsort, 'copyto': copyto, 'putmask': putmask,
    'where': where, 'nonzero': nonzero,
    'asarray': asarray, 'asanyarray': asanyarray,
    'ascontiguousarray': ascontiguousarray, 'asfortranarray': asfortranarray,
    'may_share_memory': may_share_memory, 'shares_memory': shares_memory,
    'unravel_index': unravel_index, 'ravel_multi_index': ravel_multi_index,
    'packbits': packbits, 'unpackbits': unpackbits,
    'busday_offset': busday_offset, 'busday_count': busday_count,
    'is_busday': is_busday, 'datetime_data': datetime_data,
    'datetime_as_string': datetime_as_string,
    'normalize_axis_index': normalize_axis_index,
    'scalar': scalar, '_reconstruct': _reconstruct,
    'nested_iters': nested_iters,
    'dragon4_positional': dragon4_positional,
    'dragon4_scientific': dragon4_scientific,
    'format_longfloat': format_longfloat,
    '_add_newdoc_ufunc': _add_newdoc_ufunc,
    'add_docstring': add_docstring,
    '_get_extobj_dict': _get_extobj_dict, '_make_extobj': _make_extobj,
    '_extobj_contextvar': _extobj_contextvar,
    '_get_implementing_args': _get_implementing_args,
    '_get_castingimpl': _get_castingimpl,
    '_get_sfloat_dtype': _get_sfloat_dtype,
    '_get_ndarray_c_version': _get_ndarray_c_version,
    '_populate_finfo_constants': _populate_finfo_constants,
    '_blas_supports_fpe': _blas_supports_fpe,
    '_set_madvise_hugepage': _set_madvise_hugepage,
    '_get_madvise_hugepage': _get_madvise_hugepage,
    '_set_numpy_warn_if_no_mem_policy': _set_numpy_warn_if_no_mem_policy,
    '_reload_guard': _reload_guard,
    'set_datetimeparse_function': set_datetimeparse_function,
    'set_typeDict': set_typeDict,
    'compare_chararrays': compare_chararrays,
    '_vec_string': _vec_string,
    '_discover_array_parameters': _discover_array_parameters,
    '_monotonicity': _monotonicity,
    '_place': _place,
    '_unique_hash': _unique_hash,
    '_load_from_filelike': _load_from_filelike,
    'get_handler_name': get_handler_name,
    'get_handler_version': get_handler_version,
    'from_dlpack': from_dlpack,
    'fromfile': fromfile,
    'empty_like': empty_like,
    'typeinfo': typeinfo,
    'ALLOW_THREADS': ALLOW_THREADS, 'BUFSIZE': BUFSIZE,
    'CLIP': CLIP, 'WRAP': WRAP, 'RAISE': RAISE,
    'MAXDIMS': MAXDIMS, 'MAY_SHARE_BOUNDS': MAY_SHARE_BOUNDS,
    'MAY_SHARE_EXACT': MAY_SHARE_EXACT,
    'FLOATING_POINT_SUPPORT': FLOATING_POINT_SUPPORT,
    'FPE_DIVIDEBYZERO': FPE_DIVIDEBYZERO, 'FPE_INVALID': FPE_INVALID,
    'FPE_OVERFLOW': FPE_OVERFLOW, 'FPE_UNDERFLOW': FPE_UNDERFLOW,
    'ITEM_HASOBJECT': ITEM_HASOBJECT, 'ITEM_IS_POINTER': ITEM_IS_POINTER,
    'LIST_PICKLE': LIST_PICKLE, 'NEEDS_INIT': NEEDS_INIT,
    'NEEDS_PYAPI': NEEDS_PYAPI, 'USE_GETITEM': USE_GETITEM,
    'USE_SETITEM': USE_SETITEM,
    'UFUNC_BUFSIZE_DEFAULT': UFUNC_BUFSIZE_DEFAULT,
    'UFUNC_PYVALS_NAME': UFUNC_PYVALS_NAME,
    'NAN': NAN, 'NINF': NINF, 'PINF': PINF, 'NZERO': NZERO, 'PZERO': PZERO,
    'e': e, 'pi': pi, 'euler_gamma': euler_gamma,
    '_ARRAY_API': _ARRAY_API, '_UFUNC_API': _UFUNC_API,
    '_flagdict': _flagdict, 'DATETIMEUNITS': DATETIMEUNITS,
    'tracemalloc_domain': tracemalloc_domain,
    'finfo': finfo, 'iinfo': iinfo,
}

# ===========================================================================
# Section 13 – FFT: _pocketfft_umath
# ===========================================================================
def _dft_1d(x, sign=-1):
    n = len(x)
    if n == 1:
        return list(x)
    if n & (n - 1):
        result = []
        for k in range(n):
            s = sum(x[j] * cmath.exp(sign * 2j * math.pi * k * j / n) for j in range(n))
            result.append(s)
        return result
    even = _dft_1d(x[0::2], sign)
    odd = _dft_1d(x[1::2], sign)
    T = [cmath.exp(sign * 2j * math.pi * k / n) * odd[k] for k in range(n // 2)]
    return [even[k] + T[k] for k in range(n // 2)] + \
           [even[k] - T[k] for k in range(n // 2)]

def _fft_nd(a, axes, sign=-1):
    a = array(a)
    data = [complex(x) for x in a._data]
    result = _dft_1d(data, sign)
    dt = complex128_dtype
    return _array_from_flat(result, a._shape, dt)

def _pocketfft_fft(a, n=None, axis=-1, norm=None):
    a = array(a)
    data = [complex(x) for x in a._data]
    if n is not None:
        data = data[:n] + [0j] * max(0, n - len(data))
    result = _dft_1d(data, sign=-1)
    if norm == 'ortho':
        scale = 1.0 / math.sqrt(len(result))
        result = [v * scale for v in result]
    return _array_from_flat(result, (len(result),), complex128_dtype)

def _pocketfft_ifft(a, n=None, axis=-1, norm=None):
    a = array(a)
    data = [complex(x) for x in a._data]
    if n is not None:
        data = data[:n] + [0j] * max(0, n - len(data))
    result = _dft_1d(data, sign=+1)
    m = len(result)
    if norm != 'ortho':
        result = [v / m for v in result]
    else:
        scale = 1.0 / math.sqrt(m)
        result = [v * scale for v in result]
    return _array_from_flat(result, (m,), complex128_dtype)

def _pocketfft_rfft_n_even(a, n=None, axis=-1, norm=None):
    full = _pocketfft_fft(a, n, axis, norm)
    half = (n or a.size) // 2 + 1
    return _array_from_flat(full._data[:half], (half,), complex128_dtype)

def _pocketfft_rfft_n_odd(a, n=None, axis=-1, norm=None):
    full = _pocketfft_fft(a, n, axis, norm)
    half = ((n or a.size) + 1) // 2
    return _array_from_flat(full._data[:half], (half,), complex128_dtype)

def _pocketfft_irfft(a, n=None, axis=-1, norm=None):
    a = array(a)
    data_c = a._data
    if n is None:
        n = (len(data_c) - 1) * 2
    full = list(data_c)
    while len(full) < n:
        full.append(full[n - len(full)].conjugate())
    result = _dft_1d(full, sign=+1)[:n]
    if norm != 'ortho':
        result = [v.real / n for v in result]
    else:
        result = [v.real / math.sqrt(n) for v in result]
    return _array_from_flat(result, (n,), float64_dtype)

_POCKETFFT_ATTRS = {
    'fft': ufunc('fft', 1, 1, lambda x: x),
    'ifft': ufunc('ifft', 1, 1, lambda x: x),
    'rfft_n_even': ufunc('rfft_n_even', 1, 1, lambda x: x),
    'rfft_n_odd': ufunc('rfft_n_odd', 1, 1, lambda x: x),
    'irfft': ufunc('irfft', 1, 1, lambda x: x),
    '_fft': _pocketfft_fft, '_ifft': _pocketfft_ifft,
    '_rfft_n_even': _pocketfft_rfft_n_even, '_rfft_n_odd': _pocketfft_rfft_n_odd,
    '_irfft': _pocketfft_irfft,
}

# ===========================================================================
# Section 14 – linalg
# ===========================================================================
def _linalg_solve(a, b):
    a = array(a)
    b = array(b)
    n = a._shape[0]
    aug = [list(a._data[i * n:(i + 1) * n]) + [b._data[i]] for i in range(n)]
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        if abs(pivot) < 1e-12:
            raise ValueError('Singular matrix')
        for row in range(col + 1, n):
            factor = aug[row][col] / pivot
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (aug[i][n] - sum(aug[i][j] * x[j] for j in range(i + 1, n))) / aug[i][i]
    return _array_from_flat(x, (n,), float64_dtype)

def _linalg_det(a):
    a = array(a)
    n = a._shape[0]
    mat = [list(a._data[i * n:(i + 1) * n]) for i in range(n)]
    sign = 1.0
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if max_row != col:
            mat[col], mat[max_row] = mat[max_row], mat[col]
            sign *= -1
        if abs(mat[col][col]) < 1e-14:
            return 0.0
        for row in range(col + 1, n):
            factor = mat[row][col] / mat[col][col]
            for j in range(col, n):
                mat[row][j] -= factor * mat[col][j]
    det = sign
    for i in range(n):
        det *= mat[i][i]
    return det

def _linalg_inv(a):
    a = array(a)
    n = a._shape[0]
    mat = [list(a._data[i * n:(i + 1) * n]) + [1.0 if i == j else 0.0 for j in range(n)]
           for i in range(n)]
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(mat[r][col]))
        mat[col], mat[pivot_row] = mat[pivot_row], mat[col]
        pivot = mat[col][col]
        if abs(pivot) < 1e-14:
            raise ValueError('Singular matrix')
        for j in range(2 * n):
            mat[col][j] /= pivot
        for row in range(n):
            if row != col:
                factor = mat[row][col]
                for j in range(2 * n):
                    mat[row][j] -= factor * mat[col][j]
    result = [mat[i][n + j] for i in range(n) for j in range(n)]
    return _array_from_flat(result, (n, n), float64_dtype)

def _linalg_norm(a, ord=None, axis=None, keepdims=False):
    return norm(a, ord, axis, keepdims)

def _make_linalg_stub(name):
    def _fn(*args, **kw):
        raise NotImplementedError(
            f'numpy.linalg.{name}() requires LAPACK; '
            f'use scipy or a full numpy build.')
    _fn.__name__ = name
    return _fn

_UMATH_LINALG_ATTRS = {
    '_ilp64': False,
    'solve': ufunc('solve', 2, 1, lambda a, b: 0),
    'solve1': ufunc('solve1', 2, 1, lambda a, b: 0),
    'inv': ufunc('inv', 1, 1, lambda a: 0),
    'det': ufunc('det', 1, 1, lambda a: 0),
    'slogdet': ufunc('slogdet', 1, 2, lambda a: (0, 0)),
    'eig': ufunc('eig', 1, 2, lambda a: (0, 0)),
    'eigvals': ufunc('eigvals', 1, 1, lambda a: 0),
    'eigh_lo': ufunc('eigh_lo', 1, 2, lambda a: (0, 0)),
    'eigh_up': ufunc('eigh_up', 1, 2, lambda a: (0, 0)),
    'eigvalsh_lo': ufunc('eigvalsh_lo', 1, 1, lambda a: 0),
    'eigvalsh_up': ufunc('eigvalsh_up', 1, 1, lambda a: 0),
    'cholesky_lo': ufunc('cholesky_lo', 1, 1, lambda a: 0),
    'cholesky_up': ufunc('cholesky_up', 1, 1, lambda a: 0),
    'svd': ufunc('svd', 1, 3, lambda a: (0, 0, 0)),
    'svd_f': ufunc('svd_f', 1, 3, lambda a: (0, 0, 0)),
    'svd_s': ufunc('svd_s', 1, 3, lambda a: (0, 0, 0)),
    'lstsq': ufunc('lstsq', 2, 4, lambda a, b: (0, 0, 0, 0)),
    'qr_complete': ufunc('qr_complete', 1, 2, lambda a: (0, 0)),
    'qr_reduced': ufunc('qr_reduced', 1, 2, lambda a: (0, 0)),
    'qr_r_raw': ufunc('qr_r_raw', 1, 2, lambda a: (0, 0)),
    '_solve': _linalg_solve, '_det': _linalg_det, '_inv': _linalg_inv,
    '_norm': _linalg_norm,
}

class _LapackError(Exception):
    pass

def _lapack_stub(name):
    def fn(*args, **kw):
        raise NotImplementedError(
            f'lapack_lite.{name}() requires LAPACK libraries.')
    fn.__name__ = name
    return fn

_LAPACK_LITE_ATTRS = {
    'LapackError': _LapackError, '_ilp64': False,
    'dgelsd': _lapack_stub('dgelsd'), 'dgeqrf': _lapack_stub('dgeqrf'),
    'dorgqr': _lapack_stub('dorgqr'), 'xerbla': _lapack_stub('xerbla'),
    'zgelsd': _lapack_stub('zgelsd'), 'zgeqrf': _lapack_stub('zgeqrf'),
    'zungqr': _lapack_stub('zungqr'),
}

# ===========================================================================
# Section 15 – Random module infrastructure
# ===========================================================================
class SeedlessSeedSequence:
    def generate_state(self, n, dtype=None): return [0] * n
    def spawn(self, n): return [SeedlessSeedSequence() for _ in range(n)]

SeedlessSequence = SeedlessSeedSequence

class SeedSequence:
    def __init__(self, entropy=None, *, spawn_key=(), pool_size=4):
        self._entropy = entropy
        self._rng = _random.Random(entropy)
        self.entropy = entropy
        self.spawn_key = spawn_key
        self.pool_size = pool_size
    def generate_state(self, n, dtype=None):
        return [self._rng.getrandbits(32) for _ in range(n)]
    def spawn(self, n):
        return [SeedSequence(self._rng.getrandbits(128)) for _ in range(n)]
    @property
    def state(self):
        return {'s': self._rng.getstate()}

class ISeedSequence(abc.ABC):
    @abc.abstractmethod
    def generate_state(self, n, dtype=None): ...

class ISpawnableSeedSequence(ISeedSequence):
    @abc.abstractmethod
    def spawn(self, n): ...

DECIMAL_RE = re.compile(r'\d+')

def _int_to_uint32_array(n):
    result = []
    while n:
        result.append(n & 0xFFFFFFFF)
        n >>= 32
    return result or [0]

def _coerce_to_uint32_array(x):
    if isinstance(x, int):
        return _int_to_uint32_array(x)
    return [int(v) & 0xFFFFFFFF for v in x]

class BitGenerator:
    def __init__(self, seed=None):
        if isinstance(seed, SeedSequence):
            state = seed.generate_state(4)
            seed = state[0]
        self._rng = _random.Random(seed)
        self._seed = seed
        self.lock = None
        try:
            import threading
            self.lock = threading.RLock()
        except ImportError:
            pass
    @property
    def state(self):
        return {'bit_generator': type(self).__name__,
                'state': {'rng_state': self._rng.getstate()}}
    @state.setter
    def state(self, value):
        if 'state' in value and 'rng_state' in value['state']:
            self._rng.setstate(value['state']['rng_state'])
    def random_raw(self, size=None, output=True):
        if size is None:
            return self._rng.getrandbits(64)
        n = _prod(size) if isinstance(size, tuple) else size
        data = [self._rng.getrandbits(64) for _ in range(n)]
        shape = size if isinstance(size, tuple) else (size,)
        return _array_from_flat(data, shape, uint64_dtype)
    def jumped(self, jumps=1):
        new = type(self)(self._rng.getrandbits(64))
        return new
    def spawn(self, n_children):
        return [type(self)(self._rng.getrandbits(128)) for _ in range(n_children)]
    def _generate_uint64(self, n):
        return [self._rng.getrandbits(64) for _ in range(n)]
    def _generate_double(self, n):
        return [self._rng.random() for _ in range(n)]

class MT19937(BitGenerator):
    def __init__(self, seed=None):
        super().__init__(seed)
    @property
    def state(self):
        return {'bit_generator': 'MT19937',
                'state': {'key': self._generate_uint64(624), 'pos': 624}}
    @state.setter
    def state(self, v):
        pass

class PCG64(BitGenerator):
    def __init__(self, seed=None):
        super().__init__(seed)
    @property
    def state(self):
        return {'bit_generator': 'PCG64',
                'state': {'state': self._rng.getrandbits(128), 'inc': 1}}
    @state.setter
    def state(self, v):
        pass

class PCG64DXSM(PCG64):
    pass

class Philox(BitGenerator):
    def __init__(self, seed=None, counter=None, key=None):
        super().__init__(seed)
    @property
    def state(self):
        return {'bit_generator': 'Philox',
                'state': {'counter': [0]*4, 'key': [0]*2, 'buffer': [0]*4, 'buffer_pos': 4}}
    @state.setter
    def state(self, v): pass
    @property
    def ctypes(self): return None
    @property
    def cffi(self): return None

class SFC64(BitGenerator):
    def __init__(self, seed=None):
        super().__init__(seed)

_BIT_GENERATOR_ATTRS = {
    'BitGenerator': BitGenerator, 'SeedSequence': SeedSequence,
    'ISeedSequence': ISeedSequence, 'ISpawnableSeedSequence': ISpawnableSeedSequence,
    'SeedlessSeedSequence': SeedlessSeedSequence, 'SeedlessSequence': SeedlessSequence,
    'DECIMAL_RE': DECIMAL_RE, '_int_to_uint32_array': _int_to_uint32_array,
    '_coerce_to_uint32_array': _coerce_to_uint32_array,
    'randbits': _random.getrandbits, 'cycle': itertools.cycle,
    'abc': abc, 're': re, 'np': None,
}

_MT19937_ATTRS = {'MT19937': MT19937, 'np': None, 'operator': operator}
_PCG64_ATTRS = {'PCG64': PCG64, 'PCG64DXSM': PCG64DXSM, 'np': None}
_PHILOX_ATTRS = {'Philox': Philox, 'np': None}
_SFC64_ATTRS = {'SFC64': SFC64, 'np': None}
_BOUNDED_INT_ATTRS = {'np': None}
_COMMON_ATTRS = {'interface': None, 'namedtuple': namedtuple, 'np': None, 'sys': sys}

class Generator:
    def __init__(self, bit_generator):
        self.bit_generator = bit_generator
        self._rng = bit_generator._rng
    def integers(self, low, high=None, size=None, dtype=None, endpoint=False):
        if high is None:
            low, high = 0, low
        if endpoint:
            high += 1
        if size is None:
            return self._rng.randint(int(low), int(high) - 1)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.randint(int(low), int(high) - 1) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), int64_dtype)
    def random(self, size=None, dtype=None, out=None):
        if size is None:
            return self._rng.random()
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.random() for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def standard_normal(self, size=None, dtype=None, out=None):
        if size is None:
            return self._rng.gauss(0, 1)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.gauss(0, 1) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def normal(self, loc=0.0, scale=1.0, size=None):
        if size is None:
            return self._rng.gauss(loc, scale)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.gauss(loc, scale) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def uniform(self, low=0.0, high=1.0, size=None):
        if size is None:
            return self._rng.uniform(low, high)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.uniform(low, high) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def choice(self, a, size=None, replace=True, p=None, axis=0, shuffle=True):
        a = list(array(a)._data)
        if size is None:
            return self._rng.choice(a)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        if replace:
            data = [self._rng.choice(a) for _ in range(n)]
        else:
            data = self._rng.sample(a, n)
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def shuffle(self, x, axis=0):
        x = array(x)
        self._rng.shuffle(x._data)
    def permutation(self, x):
        if isinstance(x, int):
            data = list(range(x))
            self._rng.shuffle(data)
            return _array_from_flat(data, (x,), int64_dtype)
        x = array(x).copy()
        self._rng.shuffle(x._data)
        return x
    def permuted(self, x, axis=None, out=None):
        return self.permutation(x)
    def exponential(self, scale=1.0, size=None):
        if size is None:
            return self._rng.expovariate(1.0 / scale)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.expovariate(1.0 / scale) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def poisson(self, lam=1.0, size=None):
        def _poisson(lam):
            L = math.exp(-lam)
            k, p = 0, 1.0
            while p > L:
                k += 1
                p *= self._rng.random()
            return k - 1
        if size is None:
            return _poisson(lam)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [_poisson(lam) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), int64_dtype)
    def binomial(self, n, p, size=None):
        if size is None:
            return sum(1 for _ in range(n) if self._rng.random() < p)
        sz = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [sum(1 for _ in range(n) if self._rng.random() < p) for _ in range(sz)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), int64_dtype)
    def beta(self, a, b, size=None):
        if size is None:
            return self._rng.betavariate(a, b)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.betavariate(a, b) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def gamma(self, shape_param, scale=1.0, size=None):
        if size is None:
            return self._rng.gammavariate(shape_param, scale)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [self._rng.gammavariate(shape_param, scale) for _ in range(n)]
        sh = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(sh), float64_dtype)
    def dirichlet(self, alpha, size=None):
        samples = [self._rng.gammavariate(a, 1) for a in alpha]
        s = sum(samples)
        samples = [v / s for v in samples]
        if size is None:
            return _array_from_flat(samples, (len(samples),), float64_dtype)
        return _array_from_flat(samples, (len(samples),), float64_dtype)
    def multinomial(self, n, pvals, size=None):
        cum = []
        c = 0
        for p in pvals:
            c += p
            cum.append(c)
        result = [0] * len(pvals)
        for _ in range(n):
            r = self._rng.random()
            for i, c in enumerate(cum):
                if r <= c:
                    result[i] += 1
                    break
        return _array_from_flat(result, (len(result),), int64_dtype)
    def multivariate_normal(self, mean, cov, size=None):
        k = len(mean)
        if size is None:
            data = [mean[i] + self._rng.gauss(0, math.sqrt(abs(cov[i][i])))
                    for i in range(k)]
            return _array_from_flat(data, (k,), float64_dtype)
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        all_data = []
        for _ in range(n):
            all_data.extend([mean[i] + self._rng.gauss(0, math.sqrt(abs(cov[i][i])))
                             for i in range(k)])
        shape = (size if isinstance(size, int) else _prod(size), k)
        return _array_from_flat(all_data, shape, float64_dtype)
    def standard_cauchy(self, size=None):
        if size is None:
            u = self._rng.random()
            return math.tan(math.pi * (u - 0.5))
        n = _prod(size) if isinstance(size, (tuple, list)) else size
        data = [math.tan(math.pi * (self._rng.random() - 0.5)) for _ in range(n)]
        shape = size if isinstance(size, (tuple, list)) else (size,)
        return _array_from_flat(data, tuple(shape), float64_dtype)
    def standard_exponential(self, size=None):
        return self.exponential(1.0, size)
    def standard_gamma(self, shape_param, size=None):
        return self.gamma(shape_param, 1.0, size)
    def __repr__(self):
        return f'Generator({type(self.bit_generator).__name__})'
    @property
    def __class__(self):
        return Generator

def default_rng(seed=None):
    if seed is None or isinstance(seed, (int, float)):
        bg = PCG64(int(seed) if seed is not None else None)
    elif isinstance(seed, BitGenerator):
        bg = seed
    elif isinstance(seed, SeedSequence):
        bg = PCG64(seed.generate_state(1)[0])
    else:
        bg = PCG64(int(seed))
    return Generator(bg)

_GENERATOR_ATTRS = {
    'Generator': Generator, 'PCG64': PCG64, 'default_rng': default_rng,
    'Sequence': _cabc.Sequence, 'normalize_axis_index': normalize_axis_index,
    'np': None, 'operator': operator, 'warnings': None,
}

class RandomState:
    def __init__(self, seed=None):
        self._mt19937 = MT19937(seed)
        self._rng = self._mt19937._rng
    def _gen(self):
        return Generator(self._mt19937)
    def seed(self, seed=None):
        self._mt19937 = MT19937(seed)
        self._rng = self._mt19937._rng
    def get_state(self, legacy=True):
        return {'bit_generator': 'MT19937',
                'state': {'key': list(range(624)), 'pos': 0},
                'has_gauss': 0, 'cached_gaussian': 0.0}
    def set_state(self, state): pass
    def random_sample(self, size=None): return self._gen().random(size)
    def random(self, size=None): return self._gen().random(size)
    def rand(self, *args): return self._gen().random(args if args else None)
    def randn(self, *args): return self._gen().standard_normal(args if args else None)
    def randint(self, low, high=None, size=None, dtype=None):
        return self._gen().integers(low, high, size)
    def random_integers(self, low, high=None, size=None):
        return self._gen().integers(low, (high or low) + 1, size)
    def choice(self, a, size=None, replace=True, p=None):
        return self._gen().choice(a, size, replace, p)
    def shuffle(self, x): return self._gen().shuffle(x)
    def permutation(self, x): return self._gen().permutation(x)
    def uniform(self, low=0.0, high=1.0, size=None):
        return self._gen().uniform(low, high, size)
    def normal(self, loc=0.0, scale=1.0, size=None):
        return self._gen().normal(loc, scale, size)
    def standard_normal(self, size=None): return self._gen().standard_normal(size)
    def exponential(self, scale=1.0, size=None): return self._gen().exponential(scale, size)
    def poisson(self, lam=1.0, size=None): return self._gen().poisson(lam, size)
    def binomial(self, n, p, size=None): return self._gen().binomial(n, p, size)
    def beta(self, a, b, size=None): return self._gen().beta(a, b, size)
    def gamma(self, shape, scale=1.0, size=None): return self._gen().gamma(shape, scale, size)
    def dirichlet(self, alpha, size=None): return self._gen().dirichlet(alpha, size)
    def multinomial(self, n, pvals, size=None): return self._gen().multinomial(n, pvals, size)
    def multivariate_normal(self, mean, cov, size=None):
        return self._gen().multivariate_normal(mean, cov, size)
    def bytes(self, length): return bytes([self._rng.randint(0, 255) for _ in range(length)])
    def __repr__(self): return 'RandomState(MT19937)'
    @property
    def __class__(self): return RandomState

_rand = RandomState()

def _mtrand_setter(seed=None): _rand.seed(seed)

_MTRAND_ATTRS = {
    'RandomState': RandomState, 'Sequence': _cabc.Sequence, '_MT19937': MT19937,
    'seed': _rand.seed, 'get_state': _rand.get_state, 'set_state': _rand.set_state,
    'random_sample': _rand.random, 'random': _rand.random, 'ranf': _rand.random,
    'rand': lambda *shape: _rand.random(shape if shape else None),
    'randn': lambda *shape: _rand.normal(0, 1, shape if shape else None),
    'randint': _rand.randint,
    'random_integers': lambda low, high=None, size=None: _rand.randint(
        low, (high or low) + 1, size),
    'choice': _rand.choice, 'shuffle': _rand.shuffle, 'permutation': _rand.permutation,
    'bytes': lambda n: bytes(_rand.randint(0, 256, n).tolist()),
    'normal': _rand.normal, 'uniform': _rand.uniform,
    'standard_normal': _rand.standard_normal,
    'standard_exponential': lambda size=None: _rand.exponential(1.0, size),
    'standard_gamma': lambda shape, size=None: _rand.gamma(shape, 1.0, size),
    'standard_cauchy': lambda size=None: _rand.standard_normal(size) / _rand.standard_normal(size),
    'exponential': _rand.exponential, 'poisson': _rand.poisson,
    'binomial': _rand.binomial,
    'negative_binomial': lambda n, p, size=None: zeros(size or 1),
    'geometric': lambda p, size=None: zeros(size or 1),
    'hypergeometric': lambda ngood, nbad, nsample, size=None: zeros(size or 1),
    'logseries': lambda p, size=None: zeros(size or 1),
    'multinomial': _rand.multinomial, 'multivariate_normal': _rand.multivariate_normal,
    'dirichlet': _rand.dirichlet, 'beta': _rand.beta, 'gamma': _rand.gamma,
    'f': lambda dfnum, dfden, size=None: zeros(size or 1),
    'chisquare': lambda df, size=None: zeros(size or 1),
    'noncentral_chisquare': lambda df, nonc, size=None: zeros(size or 1),
    'noncentral_f': lambda dfnum, dfden, nonc, size=None: zeros(size or 1),
    'vonmises': lambda mu, kappa, size=None: zeros(size or 1),
    'pareto': lambda a, size=None: zeros(size or 1),
    'weibull': lambda a, size=None: zeros(size or 1),
    'power': lambda a, size=None: zeros(size or 1),
    'laplace': lambda loc=0.0, scale=1.0, size=None: _rand.normal(loc, scale, size),
    'gumbel': lambda loc=0.0, scale=1.0, size=None: zeros(size or 1),
    'logistic': lambda loc=0.0, scale=1.0, size=None: zeros(size or 1),
    'lognormal': lambda mean=0.0, sigma=1.0, size=None: zeros(size or 1),
    'wald': lambda mean, scale, size=None: zeros(size or 1),
    'triangular': lambda left, mode, right, size=None: zeros(size or 1),
    'rayleigh': lambda scale=1.0, size=None: zeros(size or 1),
    'get_bit_generator': lambda: MT19937(), 'set_bit_generator': lambda bg: None,
    'zipf': lambda a, size=None: zeros(size or 1),
    'tomaxint': lambda size=None: _rand.randint(0, 2**31, size),
}

# ===========================================================================
# Section 16 – _simd stub
# ===========================================================================
class _SIMDModule:
    baseline = {}
    targets = {}
    def clear_floatstatus(self): pass
    def get_floatstatus(self): return 0

_SIMD_ATTRS = {
    'X86_V3': None, 'X86_V4': None, 'baseline': {}, 'targets': {},
    'clear_floatstatus': lambda: None, 'get_floatstatus': lambda: 0,
}

# ===========================================================================
# Section 17 – Test module stubs
# ===========================================================================
_TEST_STUB_ATTRS = {}

# ===========================================================================
# Section 18 – install() – register all shims in sys.modules
# ===========================================================================
def install(also_legacy_core=True):
    _make_module('numpy._core._multiarray_umath', **_MULTIARRAY_UMATH_ATTRS)
    _make_module('numpy._core._multiarray_tests', **_TEST_STUB_ATTRS)
    _make_module('numpy._core._operand_flag_tests', **_TEST_STUB_ATTRS)
    _make_module('numpy._core._rational_tests', **_TEST_STUB_ATTRS)
    _make_module('numpy._core._simd', **_SIMD_ATTRS)
    _make_module('numpy._core._struct_ufunc_tests', **_TEST_STUB_ATTRS)
    _make_module('numpy._core._umath_tests', **_TEST_STUB_ATTRS)
    _make_module('numpy.fft._pocketfft_umath', **_POCKETFFT_ATTRS)
    _make_module('numpy.linalg._umath_linalg', **_UMATH_LINALG_ATTRS)
    _make_module('numpy.linalg.lapack_lite', **_LAPACK_LITE_ATTRS)
    _make_module('numpy.random.bit_generator', **_BIT_GENERATOR_ATTRS)
    _make_module('numpy.random._bounded_integers', **_BOUNDED_INT_ATTRS)
    _make_module('numpy.random._common', **_COMMON_ATTRS)
    _make_module('numpy.random._generator', **_GENERATOR_ATTRS)
    _make_module('numpy.random._mt19937', **_MT19937_ATTRS)
    _make_module('numpy.random._pcg64', **_PCG64_ATTRS)
    _make_module('numpy.random._philox', **_PHILOX_ATTRS)
    _make_module('numpy.random._sfc64', **_SFC64_ATTRS)
    _make_module('numpy.random.mtrand', **_MTRAND_ATTRS)

    import types as _types
    _dtypes_mod = _types.ModuleType('numpy.dtypes')
    _dtypes_mod.__package__ = 'numpy'
    _dtypes_mod.__spec__ = None
    _dtype_map = {
        'BoolDType': bool_, 'Int8DType': int8, 'UInt8DType': uint8,
        'Int16DType': int16, 'UInt16DType': uint16, 'Int32DType': int32,
        'IntDType': int32, 'UInt32DType': uint32, 'UIntDType': uint32,
        'Int64DType': int64, 'UInt64DType': uint64, 'LongLongDType': int64,
        'ULongLongDType': uint64, 'Float16DType': float16, 'Float32DType': float32,
        'Float64DType': float64, 'LongDoubleDType': float64, 'Complex64DType': complex64,
        'Complex128DType': complex128, 'CLongDoubleDType': complex128,
        'ObjectDType': object_, 'BytesDType': bytes_, 'StrDType': str_,
        'VoidDType': void, 'DateTime64DType': datetime64, 'TimeDelta64DType': timedelta64,
    }
    for _name, _cls in _dtype_map.items():
        setattr(_dtypes_mod, _name, _cls)
    _dtypes_mod.__all__ = list(_dtype_map.keys())
    sys.modules['numpy.dtypes'] = _dtypes_mod

    if also_legacy_core:
        _make_module('numpy.core._multiarray_umath', **_MULTIARRAY_UMATH_ATTRS)
        _make_module('numpy.core._multiarray_tests', **_TEST_STUB_ATTRS)
        _make_module('numpy.core._operand_flag_tests', **_TEST_STUB_ATTRS)
        _make_module('numpy.core._rational_tests', **_TEST_STUB_ATTRS)
        _make_module('numpy.core._struct_ufunc_tests', **_TEST_STUB_ATTRS)
        _make_module('numpy.core._umath_tests', **_TEST_STUB_ATTRS)
        _make_module('numpy.core._simd', **_SIMD_ATTRS)

    # Post-install: export iinfo/finfo to numpy namespace
    _np = sys.modules.get('numpy')
    if _np is not None:
        _np.iinfo = iinfo
        _np.finfo = finfo
        _np.dtype = dtype
        _np.ndarray = ndarray
        _core = sys.modules.get('numpy._core')
        if _core is not None:
            _core.iinfo = iinfo
            _core.finfo = finfo
            _core.dtype = dtype
            _core.ndarray = ndarray

    for modname, mod in list(sys.modules.items()):
        parts = modname.split('.')
        if len(parts) > 1:
            parent_name = '.'.join(parts[:-1])
            child_name = parts[-1]
            parent = sys.modules.get(parent_name)
            if parent is not None and not hasattr(parent, child_name):
                try:
                    setattr(parent, child_name, mod)
                except (AttributeError, TypeError):
                    pass

    return {k: sys.modules[k] for k in sys.modules if 'numpy' in k and
            any(s in k for s in ('_multiarray_umath', '_pocketfft', '_umath_linalg',
                                  'lapack_lite', 'bit_generator', '_mt19937',
                                  '_pcg64', '_philox', '_sfc64', 'mtrand',
                                  '_bounded_int', '_common', '_generator', '_simd'))}

install()

class _SubmoduleAttrFixer:
    def find_module(self, name, path=None): return None
    def find_spec(self, fullname, path, target=None):
        if fullname == 'numpy' or fullname.startswith('numpy.'):
            _ensure_submodule_attrs()
        return None

sys.meta_path.insert(0, _SubmoduleAttrFixer())

__all__ = [
    'install', 'ndarray', 'dtype', 'ufunc', 'broadcast', 'flatiter', 'nditer', 'flagsobj',
    'generic', 'number', 'integer', 'signedinteger', 'unsignedinteger',
    'inexact', 'floating', 'complexfloating', 'flexible', 'character',
    'bool_', 'int8', 'int16', 'int32', 'int64', 'int_', 'intp',
    'uint8', 'uint16', 'uint32', 'uint64',
    'float16', 'float32', 'float64', 'float_',
    'complex64', 'complex128', 'complex_',
    'bytes_', 'str_', 'object_', 'void', 'datetime64', 'timedelta64',
    'array', 'zeros', 'ones', 'empty', 'full', 'arange', 'linspace', 'logspace',
    'eye', 'identity', 'diag', 'concatenate', 'stack', 'vstack', 'hstack',
    'dot', 'matmul', 'vdot', 'inner', 'outer', 'einsum',
    'where', 'nonzero', 'argwhere', 'sort', 'argsort', 'lexsort', 'unique', 'searchsorted',
    'sum_', 'cumsum', 'cumprod', 'diff',
    'BitGenerator', 'MT19937', 'PCG64', 'PCG64DXSM', 'Philox', 'SFC64',
    'Generator', 'RandomState', 'SeedSequence', 'default_rng',
    'pi', 'e', 'euler_gamma', 'NAN', 'PINF', 'NINF',
    'finfo', 'iinfo',
]

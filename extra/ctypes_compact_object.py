# Dylan Kenneth Eliot

"""

This file allows for compression of objects, their data, and objective state.



"""




from __future__ import annotations

import ctypes
import pickle
import struct
import zlib
from dataclasses import dataclass
from typing import Any


# Record tags
NULL = 0
FALSE = 1
TRUE = 2
INT = 3
FLOAT64 = 4
COMPLEX128 = 5
UTF8 = 6
BYTES = 7
LIST = 8
TUPLE = 9
DICT = 10
SET = 11
FROZENSET = 12
REF = 13
PICKLE = 14


def _uvarint(n: int) -> bytes:
    if n < 0:
        raise ValueError("uvarint requires a non-negative integer")
    out = bytearray()
    while n >= 0x80:
        out.append((n & 0x7F) | 0x80)
        n >>= 7
    out.append(n)
    return bytes(out)


def _read_uvarint(data: memoryview, pos: int) -> tuple[int, int]:
    value = 0
    shift = 0
    while True:
        if pos >= len(data):
            raise ValueError("truncated variable-length integer")
        b = data[pos]
        pos += 1
        value |= (b & 0x7F) << shift
        if not (b & 0x80):
            return value, pos
        shift += 7
        if shift > 63_000:
            raise ValueError("invalid variable-length integer")


def _zigzag_encode(n: int) -> int:
    return n * 2 if n >= 0 else (-n * 2) - 1


def _zigzag_decode(n: int) -> int:
    return -(n // 2) - 1 if n & 1 else n // 2


class CTypesArena:
    """Growable contiguous byte arena backed by a ctypes array."""

    def __init__(self, initial_capacity: int = 256) -> None:
        self.capacity = max(1, initial_capacity)
        self.length = 0
        self._buffer = (ctypes.c_ubyte * self.capacity)()

    def _reserve(self, additional: int) -> None:
        required = self.length + additional
        if required <= self.capacity:
            return
        new_capacity = self.capacity
        while new_capacity < required:
            new_capacity *= 2
        new_buffer = (ctypes.c_ubyte * new_capacity)()
        ctypes.memmove(new_buffer, self._buffer, self.length)
        self._buffer = new_buffer
        self.capacity = new_capacity

    def append_byte(self, value: int) -> None:
        self._reserve(1)
        self._buffer[self.length] = value
        self.length += 1

    def append(self, data: bytes | bytearray | memoryview) -> None:
        raw = bytes(data)
        self._reserve(len(raw))
        if raw:
            ctypes.memmove(
                ctypes.addressof(self._buffer) + self.length,
                raw,
                len(raw),
            )
            self.length += len(raw)

    def to_bytes(self) -> bytes:
        return ctypes.string_at(ctypes.addressof(self._buffer), self.length)


@dataclass(slots=True)
class PackedObject:
    payload: bytes
    compressed: bool
    raw_size: int

    def unpack(self) -> Any:
        raw = zlib.decompress(self.payload) if self.compressed else self.payload
        if len(raw) != self.raw_size:
            raise ValueError("corrupt packed object")
        return CompactDecoder(raw).decode()

    @property
    def stored_bytes(self) -> int:
        return len(self.payload)


class CompactEncoder:
    """
    Lossless adaptive encoder for Python object graphs.

    Native compact forms:
      None, bool, int, float, complex, str, bytes,
      list, tuple, dict, set, frozenset.

    Other pickle-compatible objects use a tagged pickle fallback.
    Repeated references and cycles are represented by reference IDs.
    """

    def __init__(self) -> None:
        self.out = CTypesArena()
        self.memo: dict[int, int] = {}

    @staticmethod
    def _is_reference_type(obj: Any) -> bool:
        return isinstance(
            obj,
            (str, bytes, bytearray, memoryview,
             list, tuple, dict, set, frozenset),
        ) or not isinstance(obj, (type(None), bool, int, float, complex))

    def _write_length_prefixed(self, raw: bytes) -> None:
        self.out.append(_uvarint(len(raw)))
        self.out.append(raw)

    def encode(self, obj: Any) -> bytes:
        self._write(obj)
        return self.out.to_bytes()

    def _write(self, obj: Any) -> None:
        if obj is None:
            self.out.append_byte(NULL)
            return
        if obj is False:
            self.out.append_byte(FALSE)
            return
        if obj is True:
            self.out.append_byte(TRUE)
            return

        if self._is_reference_type(obj):
            object_id = id(obj)
            if object_id in self.memo:
                self.out.append_byte(REF)
                self.out.append(_uvarint(self.memo[object_id]))
                return
            self.memo[object_id] = len(self.memo)

        if isinstance(obj, int):
            self.out.append_byte(INT)
            self.out.append(_uvarint(_zigzag_encode(obj)))
        elif isinstance(obj, float):
            self.out.append_byte(FLOAT64)
            self.out.append(struct.pack("<d", obj))
        elif isinstance(obj, complex):
            self.out.append_byte(COMPLEX128)
            self.out.append(struct.pack("<dd", obj.real, obj.imag))
        elif isinstance(obj, str):
            self.out.append_byte(UTF8)
            self._write_length_prefixed(obj.encode("utf-8"))
        elif isinstance(obj, (bytes, bytearray, memoryview)):
            self.out.append_byte(BYTES)
            self._write_length_prefixed(bytes(obj))
        elif isinstance(obj, list):
            self.out.append_byte(LIST)
            self.out.append(_uvarint(len(obj)))
            for item in obj:
                self._write(item)
        elif isinstance(obj, tuple):
            self.out.append_byte(TUPLE)
            self.out.append(_uvarint(len(obj)))
            for item in obj:
                self._write(item)
        elif isinstance(obj, dict):
            self.out.append_byte(DICT)
            self.out.append(_uvarint(len(obj)))
            for key, value in obj.items():
                self._write(key)
                self._write(value)
        elif isinstance(obj, set):
            self.out.append_byte(SET)
            self.out.append(_uvarint(len(obj)))
            for item in obj:
                self._write(item)
        elif isinstance(obj, frozenset):
            self.out.append_byte(FROZENSET)
            self.out.append(_uvarint(len(obj)))
            for item in obj:
                self._write(item)
        else:
            raw = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
            self.out.append_byte(PICKLE)
            self._write_length_prefixed(raw)


class CompactDecoder:
    def __init__(self, raw: bytes) -> None:
        self.data = memoryview(raw)
        self.pos = 0
        self.memo: list[Any] = []

    def _byte(self) -> int:
        if self.pos >= len(self.data):
            raise ValueError("truncated input")
        b = self.data[self.pos]
        self.pos += 1
        return b

    def _varint(self) -> int:
        value, self.pos = _read_uvarint(self.data, self.pos)
        return value

    def _take(self, size: int) -> bytes:
        end = self.pos + size
        if end > len(self.data):
            raise ValueError("truncated input")
        raw = bytes(self.data[self.pos:end])
        self.pos = end
        return raw

    def _blob(self) -> bytes:
        return self._take(self._varint())

    def decode(self) -> Any:
        value = self._read()
        if self.pos != len(self.data):
            raise ValueError("trailing bytes after packed object")
        return value

    def _read(self) -> Any:
        tag = self._byte()

        if tag == NULL:
            return None
        if tag == FALSE:
            return False
        if tag == TRUE:
            return True
        if tag == INT:
            return _zigzag_decode(self._varint())
        if tag == FLOAT64:
            return struct.unpack("<d", self._take(8))[0]
        if tag == COMPLEX128:
            real, imag = struct.unpack("<dd", self._take(16))
            return complex(real, imag)
        if tag == REF:
            ref = self._varint()
            try:
                return self.memo[ref]
            except IndexError as exc:
                raise ValueError("invalid reference ID") from exc

        if tag == UTF8:
            value = self._blob().decode("utf-8")
            self.memo.append(value)
            return value
        if tag == BYTES:
            value = self._blob()
            self.memo.append(value)
            return value

        if tag == LIST:
            count = self._varint()
            value: list[Any] = []
            self.memo.append(value)
            value.extend(self._read() for _ in range(count))
            return value

        if tag == DICT:
            count = self._varint()
            value: dict[Any, Any] = {}
            self.memo.append(value)
            for _ in range(count):
                key = self._read()
                item = self._read()
                value[key] = item
            return value

        if tag == SET:
            count = self._varint()
            value: set[Any] = set()
            self.memo.append(value)
            for _ in range(count):
                value.add(self._read())
            return value

        if tag == TUPLE:
            count = self._varint()
            placeholder: list[Any] = []
            memo_index = len(self.memo)
            self.memo.append(placeholder)
            value = tuple(self._read() for _ in range(count))
            self.memo[memo_index] = value
            return value

        if tag == FROZENSET:
            count = self._varint()
            placeholder: set[Any] = set()
            memo_index = len(self.memo)
            self.memo.append(placeholder)
            value = frozenset(self._read() for _ in range(count))
            self.memo[memo_index] = value
            return value

        if tag == PICKLE:
            value = pickle.loads(self._blob())
            self.memo.append(value)
            return value

        raise ValueError(f"unknown tag: {tag}")


def pack_object(obj: Any, compression_level: int = 9) -> PackedObject:
    raw = CompactEncoder().encode(obj)
    compressed = zlib.compress(raw, compression_level)

    # Compression is retained only when it actually reduces storage.
    if len(compressed) < len(raw):
        return PackedObject(compressed, True, len(raw))
    return PackedObject(raw, False, len(raw))


if __name__ == "__main__":
    shared = {"coordinates": [1, 2, 3], "label": "active"}
    original = {
        "first": shared,
        "second": shared,
        "values": [0] * 10_000,
        "mixed": (None, True, -7, 1.25, b"abc"),
    }
    original["cycle"] = original

    packed = pack_object(original)
    restored = packed.unpack()

    print("Packed bytes:", packed.stored_bytes)
    print("Compressed:", packed.compressed)
    print("Shared reference retained:",
          restored["first"] is restored["second"])
    print("Cycle retained:", restored["cycle"] is restored)

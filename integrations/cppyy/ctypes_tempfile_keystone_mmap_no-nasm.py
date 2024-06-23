# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Works for cpython-3.6 to cpython3.10 

"""
import os
import tempfile
import ctypes
from ctypes import c_int, c_void_p, CFUNCTYPE, cast, POINTER, addressof
from keystone import Ks, KS_ARCH_X86, KS_MODE_64
import cppyy

from mmap import MAP_PRIVATE, PROT_WRITE, PROT_READ, PAGESIZE, PROT_EXEC, mmap

def AssemblyFunction(assembly_code):
    ks = Ks(KS_ARCH_X86, KS_MODE_64)
    encoding, _ = ks.asm(assembly_code)
    return bytes(encoding)

c_code="""
extern "C" {
    int (*f)(int);
    int main(void) {
        return f(42);
    }
}
"""
assembly_code="""
mov eax, edi
add eax, 1
ret
"""

cppyy.cppdef(c_code)
f=tempfile.TemporaryFile()
f.write(b'\x00' * PAGESIZE)
f.seek(0)
buf = mmap(f.fileno(), PAGESIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE)
assembled_code = AssemblyFunction(assembly_code)
buf.write(assembled_code)
ctypes.CDLL(None).mprotect(ctypes.c_void_p(ctypes.addressof(ctypes.c_byte.from_buffer(buf))),PAGESIZE,PROT_READ|PROT_WRITE|PROT_EXEC)
cppyy.gbl.f=CFUNCTYPE(c_int,c_int)(cast(addressof(ctypes.c_void_p.from_buffer(buf)),c_void_p).value)
print(cppyy.gbl.main())

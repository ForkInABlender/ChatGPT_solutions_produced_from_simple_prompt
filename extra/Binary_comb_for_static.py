# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This analyzes static binaries. This can be used with tools like cheat-engine to find a segment register or function one wishes to override from the outside.
 Now looking at dockerd specifically, it has many data points. 

What it is beneficial for is when you have to know what it does without using nm, readelf, ldd, gdb or objdump.
 Because it recycles their python knowledge base, it keeps it simple for any root user to do as well as any nonroot user.
  Because it is entirely python code being used to analyze a binary, it saves time for finding the essential bits needed to use with unicorn-engine.
 It also makes it ideal to adapt to python through the raw assembly calls via ld, as, ar, or nasm.
   


"""

import os
import struct
from capstone import *


# Function to read the binary file into memory
def read_binary_file(file_path):
    with open(file_path, 'rb') as f:
        return f.read()

# Function to parse the ELF header
def parse_elf_header(binary_data):
    elf_header_format = '16sHHIQQQIHHHHHH'
    elf_header_size = struct.calcsize(elf_header_format)
    elf_header = struct.unpack(elf_header_format, binary_data[:elf_header_size])

    return {
        'e_ident': elf_header[0],
        'e_type': elf_header[1],
        'e_machine': elf_header[2],
        'e_version': elf_header[3],
        'e_entry': elf_header[4],
        'e_phoff': elf_header[5],
        'e_shoff': elf_header[6],
        'e_flags': elf_header[7],
        'e_ehsize': elf_header[8],
        'e_phentsize': elf_header[9],
        'e_phnum': elf_header[10],
        'e_shentsize': elf_header[11],
        'e_shnum': elf_header[12],
        'e_shstrndx': elf_header[13],
    }

# Function to parse section headers
def parse_section_headers(binary_data, elf_header):
    section_headers = []
    shoff = elf_header['e_shoff']
    shentsize = elf_header['e_shentsize']
    shnum = elf_header['e_shnum']

    for i in range(shnum):
        offset = shoff + i * shentsize
        section_headers.append(struct.unpack('IIQQQQIIQQ', binary_data[offset:offset + shentsize]))

    return section_headers


def get_section_names(binary_data, section_headers, shstrtab_header):
    shstrtab_offset = shstrtab_header[4]
    section_names = []

    for header in section_headers:
        name_offset = shstrtab_offset + header[0]
        name = b''
        while binary_data[name_offset] != 0:
            name += bytes([binary_data[name_offset]])
            name_offset += 1

        try:
            section_names.append(name.decode('utf-8'))
        except UnicodeDecodeError:
            section_names.append(name)

    return section_names


# Function to find a section by its type
def find_section_by_type(section_headers, section_type):
    for header in section_headers:
        if header[1] == section_type:
            return header
    return None

# Function to parse symbols from the symbol table
def parse_symbols(binary_data, symtab_header, strtab_header):
    sym_offset = symtab_header[4]
    sym_size = symtab_header[5]
    sym_entsize = symtab_header[6]
    strtab_offset = strtab_header[4]

    symbols = []

    print(f"Symbol Table Offset: {sym_offset}")
    print(f"Symbol Table Size: {sym_size}")
    print(f"Symbol Entry Size: {sym_entsize}")
    print(f"String Table Offset: {strtab_offset}")

    for i in range(0, sym_size, sym_entsize):
        entry_offset = sym_offset + i
        entry = binary_data[entry_offset:entry_offset + sym_entsize]
        
        if len(entry) != sym_entsize:
            print(f"Error: Expected symbol entry size of {sym_entsize} bytes, but got {len(entry)} bytes at offset {entry_offset}")
            continue

        try:
            st_name, st_info, st_other, st_shndx, st_value, st_size = struct.unpack('IBBHQQ', entry[:24])
        except struct.error as e:
            print(f"Unpacking error at offset {entry_offset}: {e}")
            continue

        print(f"Entry at offset {entry_offset}: st_name={st_name}, st_info={st_info}, st_other={st_other}, st_shndx={st_shndx}, st_value={st_value}, st_size={st_size}")

        # Validate st_name
        if st_name >= len(binary_data) - strtab_offset:
            print(f"Invalid st_name value: {st_name}")
            continue

        # Extract the name
        name_offset = strtab_offset + st_name
        name = b''
        while binary_data[name_offset] != 0:
            name += bytes([binary_data[name_offset]])
            name_offset += 1

        try:
            decoded_name = name.decode('utf-8')
            symbols.append(decoded_name)
        except UnicodeDecodeError as e:
            decoded_name = name #.decode('ascii')
            symbols.append(decoded_name)
            pass
            #print(f"Unicode decoding error at offset {name_offset}: {e}")
            #symbols.append(f"<invalid utf-8 name at offset {name_offset}>")

    return symbols


def find_section_by_name(section_headers, section_names, name):
    for header, section_name in zip(section_headers, section_names):
        if section_name == name:
            return header
    return None

def disassemble_section(binary_data, section_header):
    section_offset = section_header[4]
    section_size = section_header[5]
    code = binary_data[section_offset:section_offset + section_size]
    
    md = Cs(CS_ARCH_X86, CS_MODE_64)
    for i in md.disasm(code, section_header[2]):
        print("0x%x:\t%s\t%s" % (i.address, i.mnemonic, i.op_str))

# Main execution to read the binary and parse symbols
binary_data = read_binary_file('/usr/bin/dockerd')
elf_header = parse_elf_header(binary_data)
section_headers = parse_section_headers(binary_data, elf_header)
shstrtab_header = section_headers[elf_header['e_shstrndx']]
section_names = get_section_names(binary_data, section_headers, shstrtab_header)

SHT_SYMTAB = 2
SHT_STRTAB = 3

symtab_header = find_section_by_type(section_headers, SHT_SYMTAB)
strtab_header = find_section_by_type(section_headers, SHT_STRTAB)


if symtab_header and strtab_header:
    symbols = parse_symbols(binary_data, symtab_header, strtab_header)
    for symbol in symbols:
        print(symbol)
else:
    print("Symbol table or string table not found")

sections_to_find = ['.text', '.data', '.rodata']

for section_name in sections_to_find:
    section_header = find_section_by_name(section_headers, section_names, section_name)
    if section_header:
        print(f"Section {section_name}:")
        print(f"  Offset: 0x{section_header[4]:x}")
        print(f"  Size: {section_header[5]} bytes")
        if section_name == '.text':
            print(f"Disassembly of {section_name} section:")
            disassemble_section(binary_data, section_header)
    else:
        print(f"Section {section_name} not found")

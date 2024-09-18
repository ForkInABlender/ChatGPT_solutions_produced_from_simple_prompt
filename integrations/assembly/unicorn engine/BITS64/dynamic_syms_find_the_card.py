# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
locate symbols by section or binary 'linked against', for emulation purposes.

"""


from elftools.elf.elffile import ELFFile
from elftools.elf.enums import ENUM_ST_VISIBILITY

def decode_dynsym_dynstr(file_path):
    with open(file_path, 'rb') as file:
        elf = ELFFile(file)

        # Get the dynamic symbol table, dynamic string table, and version requirement sections
        dynsym = elf.get_section_by_name('.dynsym')
        dynstr = elf.get_section_by_name('.dynstr')
        verneed = elf.get_section_by_name('.gnu.version_r')

        if not dynsym or not dynstr:
            print("Dynamic symbol table or string table not found!")
            return

        # Parse version requirements if available
        dependencies = {}
        if verneed:
            for verneed_entry in verneed.iter_versions():
                library = verneed_entry[0].name  # Accessing the library name from the first element of the tuple
                for verneed_aux in verneed_entry[1]:  # Accessing the auxiliary list from the second element of the tuple
                    symbol_index = verneed_aux['vna_other']
                    dependencies[symbol_index] = library

        # Print header
        print(f"{'Index':<6} {'Name':<20} {'Address':<10} {'Size':<6} {'Type':<15} {'Bind':<10} {'Visibility':<10} {'Section':<20} {'Linked From':<30}")
        print('-' * 140)

        # Iterate over each symbol in the dynamic symbol table
        for i, symbol in enumerate(dynsym.iter_symbols()):
            name = symbol.name  # Symbol name from dynamic string table
            address = hex(symbol['st_value'])  # Symbol address
            size = symbol['st_size']  # Symbol size
            sym_type = symbol.entry.st_info.type  # Symbol type
            bind = symbol.entry.st_info.bind  # Symbol binding
            visibility_val = symbol['st_other'].visibility
            visibility = ENUM_ST_VISIBILITY.get(visibility_val, 'STV_UNKNOWN')  # Convert to visibility string using enum
            shndx = symbol['st_shndx']  # Section index

            # Determine where the symbol is linked from
            if shndx == 'SHN_UNDEF':
                linked_from = 'External'
                # Check if we know which library it is from
                if i in dependencies:
                    linked_from += f" (from {dependencies[i]})"
            else:
                linked_from = elf.get_section(shndx).name if shndx < len(list(elf.iter_sections())) else 'Unknown'

            # Print symbol details
            print(f"{i:<6} {name:<20} {address:<10} {size:<6} {sym_type:<15} {bind:<10} {visibility:<10} {shndx:<20} {linked_from:<30}")

# Replace 'your_elf_file.elf' with the path to your ELF file
decode_dynsym_dynstr('/usr/bin/docker')

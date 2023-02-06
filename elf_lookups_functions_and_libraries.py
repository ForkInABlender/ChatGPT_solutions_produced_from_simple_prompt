from elftools.elf.elffile import ELFFile
from elftools.elf.sections import SymbolTableSection

def find_sections_and_symbols(filename):
	with open(filename, 'rb') as f:
		elf = ELFFile(f)
		sections = [section.name for section in elf.iter_sections() if section.name != '']
	return sections

def find_functions_by_section_name(filename, section_name):
	with open(filename, 'rb') as f:
		elf = ELFFile(f)
		symbol_table = None
		for section in elf.iter_sections():
			if isinstance(section, SymbolTableSection):
				symbol_table = section
				break
		if symbol_table:
			section = elf.get_section_by_name(section_name)
			if section:
				print("Functions in the", section_name, "section:")
				for symbol in symbol_table.iter_symbols():
					if symbol['st_info']['type'] == 'STT_FUNC' and symbol['st_value'] >= section['sh_addr'] and symbol['st_value'] < section['sh_addr'] + section['sh_size']:
						print("  Function:", symbol.name, "0x%x" % symbol['st_value'])
						#print(symbol.name, list(map(lambda key: hex(symbol[key]) if type(symbol[key]) == type(int()) else symbol[key], list(dict(symbol.entry).keys())))[::])
						print("symbol['st_info']:", symbol['st_info'])
						print("section['sh_addr']:", "0x%x" % section['sh_addr'])
						print("section['sh_size']:", "0x%x" % section['sh_size'])
			else:
				print("Section", section_name, "not found.")
		else:
			print("Symbol table section not found.")

def find_libraries(filename):
    with open(filename, 'rb') as f:
        elf = ELFFile(f)
        for segment in elf.iter_segments():
            if segment.header['p_type'] == 'PT_DYNAMIC':
                for tag in segment.iter_tags():
                    if tag.entry.d_tag == 'DT_NEEDED':
                        print("Library:", tag.needed)

filename="/bin/bash"
print("lookup against %s:" % filename)
find_libraries(filename)
list(map(lambda section: find_functions_by_section_name(filename, section), find_sections_and_symbols(filename))).pop()

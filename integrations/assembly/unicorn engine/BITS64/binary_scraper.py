# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

How to view raw hex data of the app intended to be 'decompiled' for emulation.


"""

from elftools.elf.elffile import ELFFile

def read_elf_sections(file_path):
    with open(file_path, 'rb') as file:
        elf = ELFFile(file)
        for section in elf.iter_sections():
            print(f"Section: {section.name}")
            data = section.data()
            # Display the data in hex format, you can modify this as needed
            hex_data = '\\x'.join(f'{byte:02x}' for byte in data)
            hex_data = '\\x' + hex_data
            print(f"Data: {hex_data}...")  # Displaying first 60 characters for brevity
            print('-' * 40)

# Replace 'your_elf_file.elf' with the path to your ELF file
read_elf_sections('/bin/ls')

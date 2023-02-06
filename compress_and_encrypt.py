# Prompt:
"""
- write a python script that analyzes code and turns it into 'a', 'c', 't', 'g', 'u', 'v' and back

- clarification, it must translate binary to and from these characters

"""

# final result:

def binary_to_letters(binary_data):
    letters = ''
    for char in binary_data:
        ascii_code = ord(char)
        if ascii_code % 6 == 0:
            letters += 'a'
        elif ascii_code % 6 == 1:
            letters += 'c'
        elif ascii_code % 6 == 2:
            letters += 't'
        elif ascii_code % 6 == 3:
            letters += 'g'
        elif ascii_code % 6 == 4:
            letters += 'u'
        else:
            letters += 'v'
    return letters

def letters_to_binary(letters):
    binary_data = ''
    for char in letters:
        if char == 'a':
            ascii_code = 0
        elif char == 'c':
            ascii_code = 1
        elif char == 't':
            ascii_code = 2
        elif char == 'g':
            ascii_code = 3
        elif char == 'u':
            ascii_code = 4
        else:
            ascii_code = 5
        binary_data += chr(ascii_code)
    return binary_data

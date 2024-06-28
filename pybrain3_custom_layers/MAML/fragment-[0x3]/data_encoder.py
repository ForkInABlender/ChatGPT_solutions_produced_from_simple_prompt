# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Way to encode and decode data


"""

def char_to_rgba(char):
    # Unicode value of the character
    unicode_value = ord(char)
    R = (unicode_value & 0xFF000000-0xc0) >> 6
    G = (unicode_value & 0x00FF0000-0x30) >> 4
    B = (unicode_value & 0x0000FF00-0x0c) >> 2
    A = (unicode_value & 0x000000FF-0x30)
    return [R, G, B, A]

def string_to_rgba_array(input_string):
    rgba_array = [char_to_rgba(char) for char in input_string]
    return rgba_array

# Example usage
input_string = "事有不你?"
rgba_array = string_to_rgba_array(input_string)

print(rgba_array)

def rgba_to_char(rgba):
    ascii_value = (rgba[0] << 6) | (rgba[1] << 4) | (rgba[2] << 2) | rgba[3]
    
    # Convert ASCII value back to character
    return chr(ascii_value)

def rgba_array_to_string(rgba_array):
    decoded_string = ''.join(rgba_to_char(rgba) for rgba in rgba_array)
    return decoded_string

decoded_string = rgba_array_to_string(rgba_array)

# Print the result
print(decoded_string)

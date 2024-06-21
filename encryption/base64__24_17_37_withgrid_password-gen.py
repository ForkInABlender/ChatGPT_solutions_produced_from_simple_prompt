# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This version keeps the old passwords from being decoded. This format encrypts by those rules.
"""

import base64
import random
import numpy as np

# Define the 10x10 grid for the grid shift method

# List of old cell phone numbers
old_phone_numbers = [
    "1234567890",
    "0987654321",
    "5555555555",
    "1112223333"
]


grid = np.array([
    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],
    ['K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'],
    ['U', 'V', 'W', 'X', 'Y', 'Z', '0', '1', '2', '3'],
    ['4', '5', '6', '7', '8', '9', '!', '@', '#', '$'],
    ['%', '^', '&', '*', '(', ')', '-', '=', '+', '['],
    [']', '{', '}', '|', ':', ';', '<', '>', ',', '.'],
    ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'],
    ['k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't'],
    ['u', 'v', 'w', 'x', 'y', 'z', '?', '/', ' ', '`'],
    ['~', '_', '"', "'", '>', '<', '€', '£', '¥', '•']
])

def find_position(char):
    """Find the position of a character in the grid."""
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if grid[i, j] == char:
                return (i, j)
    return None

def shift_position(position, shift_value):
    """Calculate the new position after shifting."""
    row, col = position
    new_row = (row + shift_value) % grid.shape[0]
    new_col = (col + shift_value) % grid.shape[1]
    return (new_row, new_col)

def encrypt_message(message, initial_shift):
    """Encrypt a message using the ShiftGrid method."""
    encrypted_message = ""
    shift_value = initial_shift
    for char in message:
        pos = find_position(char)
        if pos:
            new_pos = shift_position(pos, shift_value)
            encrypted_message += grid[new_pos]
            shift_value += 1
        else:
            encrypted_message += char  # Keep characters not in grid unchanged
    return encrypted_message

# Function to encode phone numbers to base64
def encode_base64(phone_number):
    return base64.b64encode(phone_number.encode()).decode()

# Function to generate a password stub from base64 substrings
def generate_password_stub(encoded_numbers):
    def get_random_substring(encoded, length):
        if length > len(encoded):
            length = len(encoded)  # Adjust the length if it exceeds the encoded string length
        start = random.randint(0, len(encoded) - length)
        return encoded[start:start + length]
    
    password_stub = (
        get_random_substring(encoded_numbers, 24) +
        get_random_substring(encoded_numbers, 18) +
        get_random_substring(encoded_numbers, 37)
    )
    
    return password_stub

# Encode the phone numbers

encoded_numbers = [encode_base64(number) for number in old_phone_numbers]

# Generate password stubs by combining 3 randomly selected base64 substrings
def generate_password(stub_count=3, initial_shift=5):
    selected_stubs = random.sample(encoded_numbers, stub_count)
    combined_stub = ''.join(generate_password_stub(stub) for stub in selected_stubs)
    encrypted_password = encrypt_message(combined_stub, initial_shift)
    return encrypted_password

# Generate and print a password
password = generate_password()
print("Generated Password:", password)

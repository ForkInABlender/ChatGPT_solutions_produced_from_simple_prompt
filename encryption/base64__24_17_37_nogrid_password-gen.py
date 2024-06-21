# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This version keeps the old passwords from being decoded. This format encrypts by those rules.
"""

import base64
import random

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

# List of old cell phone numbers
old_phone_numbers = [
    "1234567890",
    "0987654321",
    "5555555555",
    "1112223333"
]

# Encode the phone numbers
encoded_numbers = [encode_base64(number) for number in old_phone_numbers]

# Generate password stubs by combining 3 randomly selected base64 substrings
def generate_password(stub_count=3):
    selected_stubs = random.sample(encoded_numbers, stub_count)
    combined_stub = ''.join(generate_password_stub(stub) for stub in selected_stubs)
    return combined_stub

# Generate and print a password
password = generate_password()
print("Generated Password:", password)

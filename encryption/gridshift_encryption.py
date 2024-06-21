# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""


GPT made this one per my request.

If you change the shift number, or the number of shuffles applied to the input.







"""

import numpy as np

# Define the 10x10 grid
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

def decrypt_message(message, initial_shift):
    """Decrypt a message using the ShiftGrid method."""
    decrypted_message = ""
    shift_value = initial_shift
    for char in message:
        pos = find_position(char)
        if pos:
            # Reverse the shift direction for decryption
            new_pos = shift_position(pos, -shift_value)
            decrypted_message += grid[new_pos]
            shift_value += 1
        else:
            decrypted_message += char  # Keep characters not in grid unchanged
    return decrypted_message

# Example usage
initial_shift = 5  # Example initial shift value
original_message = "Hello World!"
encrypted_message = encrypt_message(original_message, initial_shift)
decrypted_message = decrypt_message(encrypted_message, initial_shift)

print("Original Message:", original_message)
print("Encrypted Message:", encrypted_message)
print("Decrypted Message:", decrypted_message)

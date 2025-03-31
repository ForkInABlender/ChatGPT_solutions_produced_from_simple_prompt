# Dylan Kenneth Eliot

"""
Encode and decode 99.9981% of input by character and their BPE don't mean jack.

Here's the basics to raw token translation from their BPE for open source development.


This hereby makes BPE by function fully broken down from latin-1 with base64 and iso9660 triple layering; While smart, linear drawing and playing with etcha-sketches
 ain't cutting budget. 

So instead, their BPE gets decoded and translated by reverse engineering convention like-a-dece.


And now all their data is useful for reverse engineering off gpt-2, gpt-3, and some polish away from added in human error, aka information poison plus 
 negatives gotten from excessive masking and fuzzing.
"""

import base64
import numpy as np

def decode_string_to_float32(input_string):
  encoded_bytes = np.array(input_string.encode('latin-1'), dtype=np.str_)
  float_array = np.frombuffer(encoded_bytes, dtype=np.float32)
  return float_array

def encode_float32_to_string(float_array):
  byte_array = float_array.tobytes()
  raw_string = byte_array.decode('latin-1')
  return raw_string


file_path = "./9b5ad71b2ce5302211f9c61530b329a4922fc6a4"

with open(file_path, "rb") as file:
    for line in file:
        parts = line.split()  # Split by the first word
        
        print(str(int(parts[1]))+":"+str(f"{base64.b64decode(parts[0]).decode('latin-1')}".encode('latin-1'))[1:], decode_string_to_float32(str(f"{base64.b64decode(parts[0]).decode('latin-1')}".encode('latin-1'))[2:-1]))
        print(str(int(parts[1]))+":"+str(f"{base64.b64decode(parts[0]).decode('latin-1')}".encode('latin-1'))[1:], encode_float32_to_string(decode_string_to_float32(str(f"{base64.b64decode(parts[0]).decode('latin-1')}".encode('latin-1'))[2:-1])))

# Dylan Kenneth Eliot

"""
Given the created folder, 
"""

import numpy as np
import os

# Define the directory where the .npy files are stored
params_directory = "./layer_parameters"

# Function to load and display the parameters of each layer
def load_and_display_params(directory):
    # Check if the specified directory exists
    if not os.path.exists(directory):
        print(f"Directory '{directory}' does not exist.")
        return
    
    # Iterate through all .npy files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".npy"):
            # Construct the full path to the .npy file
            filepath = os.path.join(directory, filename)
            
            # Load the numpy array from the file
            params = np.load(filepath, allow_pickle=True)  # allow_pickle to handle complex objects

            # Check if params contain bytes-like objects and attempt conversion
            try:
                # If the data type is not numeric, we will convert each element appropriately
                if isinstance(params, np.ndarray):
                    if params.dtype == 'O':  # dtype 'O' means it's an object array, possibly containing non-numeric data
                        converted_params = []
                        for elem in params:
                            if isinstance(elem, (bytes, bytearray)):
                                # Convert bytes to integer
                                converted_params.append(int.from_bytes(elem, byteorder='little', signed=False))
                            elif isinstance(elem, (list, np.ndarray)):
                                # If it's a nested list or array, recursively flatten and convert
                                flattened = flatten_nested_array(elem)
                                converted_params.extend(flattened)
                            else:
                                # Directly append if it's already a number
                                converted_params.append(elem)
                        
                        # Convert the list to a numpy array
                        params = np.array(converted_params, dtype=np.float64)
                    else:
                        # Directly attempt to convert params to float64 if they're already numerical
                        params = params.astype(np.float64)
                    
                # Print the name of the layer and its numerical data
                print(f"Layer: {filename}")
                print(f"Shape: {params.shape}")
                print(f"Data:\n{list(params)}\n")
            
            except Exception as e:
                print(f"Failed to convert {filename} to float64: {e}")

# Helper function to flatten nested sequences
def flatten_nested_array(arr):
    flattened = []
    for item in arr:
        if isinstance(item, (list, np.ndarray)):
            flattened.extend(flatten_nested_array(item))
        elif isinstance(item, (bytes, bytearray)):
            # Convert bytes to integer value
            flattened.append(int.from_bytes(item, byteorder='little', signed=False))
        else:
            flattened.append(item)
    return flattened

# Call the function to load and display parameters from the specified directory
load_and_display_params(params_directory)

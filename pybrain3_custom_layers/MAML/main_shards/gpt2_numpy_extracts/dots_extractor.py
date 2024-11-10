# Dylan Kenneth Eliot

"""
How to extract the dots from GPT models


This is part 2.
"""

import os
import numpy as np
from tensorflow.python.tools import inspect_checkpoint as chkp
import tensorflow as tf

# Path to the checkpoint
checkpoint_path = "model.ckpt-362000"

# Path to the file containing layer parameters
layer_params_file = "layer_params_type_capture.txt"

# Directory to save the extracted layer parameters
output_directory = "layer_parameters"

# Create output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def extract_and_save_parameters_for_layer(layer_name):
    """
    Extract parameters for the specified layer from the checkpoint and save to a file.
    """
    print(f"\nExtracting parameters for layer: {layer_name}\n")
    try:
        tensor_value = tf.train.load_variable(checkpoint_path, layer_name)
        output_path = os.path.join(output_directory, f"{layer_name.replace('/', '_')}.py_txt")
        filen_ = open(output_path, 'w')
        filen_.write(str(list(list(a) for a in tensor_value)))
        filen_.close()
        print(tensor_value)
        #np.save(output_path, tensor_value)
        print(f"Saved {layer_name} to {output_path}")
        del tensor_value
    except Exception as e:
        print(f"Could not extract layer {layer_name}: {e}")

def get_layer_names_from_file(file_path):
    """
    Reads the layer names from the specified file.
    """
    layer_names = []
    with open(file_path, "r") as file:
        for line in file:
            # Extract layer names without the additional parameter details
            if '(' in line:  # Ensure it has relevant information
                layer_name = line.split(' ')[0].strip()
                if layer_name and '/' in layer_name:  # Check for valid layer format
                    layer_names.append(layer_name)
    return layer_names

# Get the layer names from the file
layer_names = get_layer_names_from_file(layer_params_file)

# Iterate over each layer and extract parameters
for layer_name in layer_names:
    extract_and_save_parameters_for_layer(layer_name)

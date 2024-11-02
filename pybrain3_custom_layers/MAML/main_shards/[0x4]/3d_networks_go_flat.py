# Dylan Kenneth Eliot

"""
With the right reshape to GPT2/3/3.5/4 paramaeter data, it will work even for spatial and visual processing.

Due to it handling 3d & 2d to 1d evaluation, it will be able to process 3d data such as MRI data.

"""

import numpy as np
from pybrain3.structure import Network, Module, LinearLayer, SigmoidLayer, FullConnection

def to_1d_from_3d(array_3d):
    """Flatten a 3D array to 1D."""
    return np.array([elem for layer in array_3d for row in layer for elem in row])

def to_1d_from_2d(array_2d):
    """Flatten a 2D array to 1D."""
    return np.array([elem for row in array_2d for elem in row])

class FlattenLayer(Module):
    def __init__(self, input_dim):
        """
        FlattenLayer to transform 2D or 3D input arrays into 1D arrays.
        
        :param input_dim: tuple, dimension of input (e.g., (3, 4, 5) for 3D or (4, 5) for 2D)
        """
        self.input_dim = input_dim
        output_dim = np.prod(input_dim)  # Calculate the flattened output dimension
        super().__init__(name="FlattenLayer", indim=output_dim, outdim=output_dim)

    def _forwardImplementation(self, inbuf, outbuf):
        """
        Forward pass that flattens input data.
        
        :param inbuf: np.array, input data buffer, either 2D or 3D array
        :param outbuf: np.array, output data buffer, will store flattened data
        """
        reshaped_input = np.array(inbuf).reshape(self.input_dim)  # Ensure input is in correct shape
        if len(self.input_dim) == 3:
            flattened_data = to_1d_from_3d(reshaped_input)
        elif len(self.input_dim) == 2:
            flattened_data = to_1d_from_2d(reshaped_input)
        else:
            raise ValueError("Only 2D or 3D inputs are supported.")
        
        outbuf[:] = flattened_data  # Copy flattened data to output buffer

    def _backwardImplementation(self, outerr, inerr):
        """
        Backward pass: reshapes the error signal to match input dimensions.
        
        :param outerr: Output error signal
        :param inerr: Input error buffer
        """
        inerr[:] = outerr.reshape(self.input_dim)  # Reshape output error to match input

class CustomNetwork(Network):
    def __init__(self, input_dim, hidden_size, output_size):
        super().__init__(name="CustomNetwork")
        
        # Initialize layers
        self.flatten_layer = FlattenLayer(input_dim)       # Flatten input layer
        self.hidden_layer = SigmoidLayer(hidden_size)      # Hidden layer
        self.output_layer = LinearLayer(output_size)       # Output layer
        
        # Add layers to the network
        self.addInputModule(self.flatten_layer)
        self.addModule(self.hidden_layer)
        self.addOutputModule(self.output_layer)
        
        # Define and add connections
        self.flatten_to_hidden = FullConnection(self.flatten_layer, self.hidden_layer)
        self.hidden_to_output = FullConnection(self.hidden_layer, self.output_layer)
        
        self.addConnection(self.flatten_to_hidden)
        self.addConnection(self.hidden_to_output)
        
        # Finalize the network structure
        self.sortModules()

    def activate(self, input_data):
        """
        Custom activation method to handle multi-dimensional input.
        
        :param input_data: np.array, multi-dimensional input array (e.g., 3D for images)
        :return: np.array, output of the network after passing through layers
        """
        # Forward the input through the network manually
        flattened_output = np.zeros(self.flatten_layer.outdim)
        self.flatten_layer._forwardImplementation(input_data, flattened_output)
        
        # Pass through hidden layer
        hidden_output = self.hidden_layer.activate(flattened_output)
        
        # Pass through output layer
        final_output = self.output_layer.activate(hidden_output)
        return final_output
"""
# Usage of CustomNetwork
input_dim = (40, 32)  # Define a 3D input shape for testing
hidden_size = 1280       # Define the number of hidden units
output_size = 1280       # Define the size of the output layer

# Instantiate CustomNetwork
network = CustomNetwork(input_dim, hidden_size, output_size)

# Generate a random 3D input to test the network
test_input = np.random.rand(*input_dim)

# Pass the multi-dimensional input directly to the custom network
output = network.activate(test_input)

print("Output from the custom network:", output)

"""

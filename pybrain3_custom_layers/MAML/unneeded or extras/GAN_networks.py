# Dylan Kenneth Eliot

"""
Simple GAN setup; extra and will be unused in actual AI development for imitation of the brain.

Make use as needed, though :)


Do note that at present, pytorch, tensorflow, and transformer dependent models are out of date as of today. Code using similar levels rounding will also fall into this same folder as this file.

At no future date will pybrain3 be extended to use them by this AI developer, as all components rawly and efficiently run on the CPU with plenty of memory for x86_64 & armhf/el-hf+ processors.

If 3d processing is needed, I'd use the opengl/GPU emulator and collect the data per frame as needed in realtime. 

"""

import numpy as np
from pybrain3.structure import Network, LinearLayer, SigmoidLayer, FullConnection, BiasUnit
from pybrain3.structure.modules import TanhLayer

# Refining GAN to use Custom 3D Networks

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

# Define the Generator Network
def build_generator():
    input_dim = (10, 10, 10)  # Example 3D input dimensions for the generator
    hidden_size = 512
    output_size = 28 * 28

    generator = CustomNetwork(input_dim, hidden_size, output_size, hidden_activation='tanh')
    return generator

# Define the Discriminator Network
def build_discriminator():
    input_dim = (28, 28)  # Example 2D input dimensions for real/fake images
    hidden_size = 512
    output_size = 1

    discriminator = CustomNetwork(input_dim, hidden_size, output_size, hidden_activation='sigmoid')
    return discriminator

# Train the GAN
def train_gan(generator, discriminator, epochs=10000, batch_size=128):
    for epoch in range(epochs):
        # Generate noise and fake images
        noise = np.random.normal(0, 1, (batch_size, *generator.flatten_layer.input_dim))
        fake_images = np.array([generator.activate(noise[i]) for i in range(batch_size)])

        # Generate random real images (replace this with actual data later)
        real_images = np.random.rand(batch_size, 28, 28)

        # Labels for discriminator
        real_labels = np.ones((batch_size, 1))
        fake_labels = np.zeros((batch_size, 1))

        # Train the discriminator
        d_loss_real = []
        d_loss_fake = []

        for i in range(batch_size):
            # Train on real images
            discriminator.reset()
            output_real = discriminator.activate(real_images[i])
            loss_real = np.mean((real_labels[i] - output_real) ** 2)
            d_loss_real.append(loss_real)

            # Train on fake images
            discriminator.reset()
            output_fake = discriminator.activate(fake_images[i])
            loss_fake = np.mean((fake_labels[i] - output_fake) ** 2)
            d_loss_fake.append(loss_fake)

        d_loss = 0.5 * (np.mean(d_loss_real) + np.mean(d_loss_fake))

        # Train the generator
        g_loss = []
        for i in range(batch_size):
            generator.reset()
            generated_noise = np.random.normal(0, 1, generator.flatten_layer.input_dim)
            generated_image = generator.activate(generated_noise)
            discriminator.reset()
            fake_output = discriminator.activate(generated_image)
            loss_gen = np.mean((real_labels[i] - fake_output) ** 2)
            g_loss.append(loss_gen)

        g_loss = np.mean(g_loss)

        # Print progress
        if epoch % 1000 == 0:
            print(f"Epoch {epoch}, Discriminator Loss: {d_loss}, Generator Loss: {g_loss}")

# Instantiate and train the GAN
generator = build_generator()
discriminator = build_discriminator()
train_gan(generator, discriminator)

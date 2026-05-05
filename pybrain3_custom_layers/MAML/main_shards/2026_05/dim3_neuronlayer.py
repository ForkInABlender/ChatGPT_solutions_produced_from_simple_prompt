# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition) & Google-Gemini

from pybrain3.structure.modules.neuronlayer import NeuronLayer
import numpy as np

class Dim3NeuronLayer(NeuronLayer):
    def __init__(self, in_dim, out_dim, embed_dropout, attn_dropout, activation_function, lr, weight_decay, gradient_clipping):
        super(Dim3NeuronLayer, self).__init__(in_dim, out_dim)
        self.embed_dropout = embed_dropout
        self.attn_dropout = attn_dropout
        self.name = "Dim3NeuronLayer"
        self.activation_function = activation_function
        self.lr = lr
        self.weight_decay = weight_decay
        self.gradient_clipping = gradient_clipping

        # Initialize weights and bias
        self.weights = np.random.randn(in_dim, out_dim) * 0.1
        self.bias = np.zeros(out_dim)

    def _forwardImplementation(self, inbuf, outbuf):
        # Forward pass implementation to handle 1D, 2D, and 3D inputs
        input_data = np.array(inbuf)
        original_shape = input_data.shape

        # Reshape to 3D if needed
        if input_data.ndim == 1:
            input_data = input_data.reshape(1, 1, -1)
        elif input_data.ndim == 2:
            input_data = input_data.reshape(1, *input_data.shape)

        # Process as 3D
        output_data = self._process3D(input_data)

        # Reshape back to original shape if needed
        if len(original_shape) == 1:
            output_data = output_data.flatten()
        elif len(original_shape) == 2:
            output_data = output_data.reshape(original_shape[0], -1)

        # Ensure the output buffer size matches the output data
        outbuf[:output_data.shape[0]] = output_data[:outbuf.shape[0]]

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        # Backward pass implementation to handle 1D, 2D, and 3D inputs
        out_data = np.array(outerr)
        original_shape = out_data.shape

        # Reshape to 3D if needed
        if out_data.ndim == 1:
            out_data = out_data.reshape(1, 1, -1)
        elif out_data.ndim == 2:
            out_data = out_data.reshape(1, *out_data.shape)

        # Process as 3D
        grad_input = self._backward3D(out_data)

        # Reshape back to original shape if needed
        if len(original_shape) == 1:
            grad_input = grad_input.flatten()
        elif len(original_shape) == 2:
            grad_input = grad_input.reshape(original_shape[0], -1)

        inerr[:grad_input.shape[0]] = grad_input[:inerr.shape[0]]

    def _process3D(self, input_data):
        # Processing logic for 3D input
        output = np.tensordot(input_data, self.weights, axes=([-1], [0])) + self.bias
        return self._apply_activation(output)

    def _backward3D(self, outerr):
        # Gradient calculation for 3D input
        grad_output = outerr * self._activation_derivative(outerr)
        return grad_output

    def _apply_activation(self, input_data):
        # Apply specified activation function
        if self.activation_function == 'gelu':
            return np.where(
                input_data > 20, 1.0,  # Prevent overflow for large positive values
                np.where(
                    input_data < -20, 0.0,  # Prevent overflow for large negative values
                    input_data * 0.5 * (1.0 + np.tanh(np.sqrt(2.0 / np.pi) * (input_data + 0.044715 * np.power(input_data, 3))))
                )
            )
        else:
            raise NotImplementedError("Only GELU activation is currently implemented")

    def _activation_derivative(self, input_data):
        # Derivative of activation function
        if self.activation_function == 'gelu':
            tanh_term = np.tanh(np.sqrt(2.0 / np.pi) * (input_data + 0.044715 * np.power(input_data, 3)))
            return np.where(
                input_data > 20, 0.0,  # Prevent overflow for large positive values
                np.where(
                    input_data < -20, 0.0,  # Prevent overflow for large negative values
                    0.5 * (1.0 + tanh_term) + 0.5 * input_data * (1 - tanh_term ** 2) * (np.sqrt(2.0 / np.pi) * (1 + 3 * 0.044715 * np.power(input_data, 2)))
                )
            )
        else:
            raise NotImplementedError("Only GELU derivative is currently implemented")

class Dim3InputLayer(Dim3NeuronLayer):
    def __init__(self, in_dim, embed_dropout, attn_dropout, activation_function, lr, weight_decay, gradient_clipping):
        # Use the base class constructor with out_dim set to in_dim for an input layer
        super(Dim3InputLayer, self).__init__(in_dim, in_dim, embed_dropout, attn_dropout, activation_function='gelu', lr=0, weight_decay=0, gradient_clipping=0)
        self.name = "Dim3InputLayer"

    def _forwardImplementation(self, inbuf, outbuf):
        # As an input layer, the forward implementation simply passes the data through
        input_data = np.array(inbuf)
        original_shape = input_data.shape

        # Reshape to 3D if needed
        if input_data.ndim == 1:
            input_data = input_data.reshape(1, 1, -1)
        elif input_data.ndim == 2:
            input_data = input_data.reshape(1, *input_data.shape)

        # No processing is needed for an input layer, simply copy input to output
        output_data = input_data

        # Reshape back to original shape if needed
        if len(original_shape) == 1:
            output_data = output_data.flatten()
        elif len(original_shape) == 2:
            output_data = output_data.reshape(original_shape[0], -1)

        outbuf[:output_data.shape[0]] = output_data[:outbuf.shape[0]]

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        # Backward pass implementation for input layer
        out_data = np.array(outerr)
        original_shape = out_data.shape

        # Reshape to 3D if needed
        if out_data.ndim == 1:
            out_data = out_data.reshape(1, 1, -1)
        elif out_data.ndim == 2:
            out_data = out_data.reshape(1, *out_data.shape)

        # Since this is an input layer, we simply propagate the gradients back without modification
        grad_input = out_data

        # Reshape back to original shape if needed
        if len(original_shape) == 1:
            grad_input = grad_input.flatten()
        elif len(original_shape) == 2:
            grad_input = grad_input.reshape(original_shape[0], -1)

        inerr[:grad_input.shape[0]] = grad_input[:inerr.shape[0]]

# Dylan Kenneth Eliot

""""
PRODUCT VOID: LINE ITEM RECALL ON CONTEXT -- DO NOT EMBED INTO YOUR MODEL AS THIS IS DOES NOT CURRENTLY WORK WITH DIM3LAYER DUE TO SHAPES OF OBJECTIONAL SPACE

ChatGPT by openai botched on this one. Do not implement or use until product release of this file is labelled safe for adaptation. 


ERROR ||: -- code will run before training. Attempt to train with this class breaks runtime of model. 
   STAMPED -:||: 2024/15/10 @ 2:28:15 pm 

All other scripts outside of this one for AI development are safe for public use as all code is tested by hand and then retested. All tests done by hand are for quality assurance.

WARNING { IMPLEMENT AT YOUR OWN DISCRESSION. ALL DAMAGES YOU DO BY NOT TAKING CAREFUL ACTION IN TESTING CAN LEAD TO DESTRUCTION OF ANOTHER'S WELL BUILT PRODUCT. PLEASE USE WISELY. }

Thank you for your cooperation.

""""



"""
This defines more granular control for network definition at the layer level. This gives a more granular set of controls for when one has to tune a custom FFT block.

Because the MAML I am building with it will only use 1 of these layers with `Dim3NeuronLayer` layers, each transformer block will consist of 2 hidden layers, 1 input layer, 1 output layer.

This, as an attention head layer, simplifies and reduces the dot-product space to search within any given network built on pybrain(/3) layers.


If you're looking to make use of classes like this one based on pybrain3, use `de3343/ai_mods_py3.10:moderngl_plus_pyvirtualdesktop` docker container with "docker run -it" while passing in the files from the directory
 you aim to run it from.

"""

from pybrain3.structure.modules.neuronlayer import NeuronLayer

class AttentionHead_Layer(NeuronLayer):
    def __init__(self, inputdim, outputdim  embed_dropout, attn_dropout, activation_function, lr, weight_decay, gradient_clipping, **kwargs):

        super(CustomNeuronLayer, self).__init__(outputdim, **kwargs)

        # Store additional configurations
        self.inputdim = inputdim
        self.embed_dropout = embed_dropout
        self.attn_dropout = attn_dropout
        self.activation_function = activation_function
        self.lr = lr
        self.weight_decay = weight_decay
        self.gradient_clipping = gradient_clipping

        # Initialize weights and biases
        self.weights = np.random.randn(inputdim, outputdim) * np.sqrt(2. / inputdim)
        self.biases = np.zeros(outputdim)

    def _forwardImplementation(self, inbuf, outbuf):
        # Forward pass implementation
        linear_output = np.dot(inbuf, self.weights) + self.biases
        
        # Apply activation function
        activated_output = self._activate(linear_output)

        # Apply embedding dropout if specified
        if self.embed_dropout > 0:
            dropout_mask = np.random.rand(*activated_output.shape) > self.embed_dropout
            activated_output *= dropout_mask / (1 - self.embed_dropout)

        # Apply attention dropout if specified
        if self.attn_dropout > 0:
            attn_mask = np.random.rand(*activated_output.shape) > self.attn_dropout
            activated_output *= attn_mask / (1 - self.attn_dropout)
        
        outbuf[:] = activated_output # Oops. Was returned instead of proper assignment.

    def _backwardImplementation(self, inbuf, outbuf, deltas):
        # Backward pass implementation
        gradients = np.dot(inbuf.T, deltas) / inbuf.shape[0]  # Compute the gradient of weights
        
        # Update weights and biases with learning rate and weight decay
        self.weights -= self.lr * (gradients + self.weight_decay * self.weights)
        self.biases -= self.lr * np.sum(deltas, axis=0) / inbuf.shape[0]  # Update biases
        
        # Apply gradient clipping if specified
        if self.gradient_clipping > 0:
            np.clip(self.weights, -self.gradient_clipping, self.gradient_clipping, out=self.weights)

        return np.dot(deltas, self.weights.T)

    def _activate(self, x):
        # Define activation function (e.g., GELU)
        if self.activation_function == 'gelu':
            return 0.5 * x * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * np.power(x, 3))))
        elif self.activation_function == 'relu':
            return np.maximum(0, x)  # ReLU activation
        elif self.activation_function == 'sigmoid':
            return 1 / (1 + np.exp(-x))  # Sigmoid activation
        elif self.activation_function == 'tanh':
            return np.tanh(x)  # Tanh activation
        else:
            raise ValueError(f"Activation function '{self.activation_function}' is not supported.")

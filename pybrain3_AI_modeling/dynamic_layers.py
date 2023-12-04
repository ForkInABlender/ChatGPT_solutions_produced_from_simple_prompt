# Dylan Kenneth Eliot & GPT-4-plugins (Beta Edition) after significant beating on it with a hammer..... and no beer.... ......

"""

This is a dynamically resizing neural network layer. It has an ideal layer size, a minimum input, and a maximum input. And it adjust for each ".activate"
 input array size within the caps. This is good for if you're using something like tiktoken and numpy to bridge the gap before using brian2 & rdkit with this.

The goal is to pair this with a few modified classes for AI neurons to match the function of human neurons. Of course, that would also require some clever
 even about the assembly away from token parsing only on graphics cards :) The reason to do it this was is that it makes remodeling a few networks easier.
 On top of being easier to adapt to things like tiktoken, it should be now possible to inside of docker containers model ones own AI and adapt it to flask.

https://hub.docker.com/r/de3343/ai_mods_py3.10/

The ideal container to use is "de3343/ai_mods_py3.10:tiktoken" for modeling at the moment. 

""" 



from pybrain3.structure.modules.neuronlayer import NeuronLayer
import numpy as np

class DynamicLayer(NeuronLayer):
    def __init__(self, layer_size, min_size, max_size):
        if not (min_size <= layer_size <= max_size):
            raise ValueError(f"Layer size must be between {min_size} and {max_size}")
        super(CustomNeuronLayer, self).__init__(max_size)
        self.layer_size = layer_size
        self.min_size = min_size
        self.max_size = max_size
        self.weights = np.random.randn(self.max_size, self.max_size)
        self.biases = np.random.randn(self.max_size)

    def _adjust_weights_biases(self, input_size):
        """ Adjust the weights and biases based on the input size """
        if not (self.min_size <= input_size <= self.max_size):
            raise ValueError(f"Input size must be between {self.min_size} and {self.max_size}")
        adjusted_weights = self.weights[:input_size, :input_size]
        adjusted_biases = self.biases[:input_size]
        return adjusted_weights, adjusted_biases

    def _forwardImplementation(self, inbuf, outbuf):
        adjusted_weights, adjusted_biases = self._adjust_weights_biases(len(inbuf))
        outbuf[:] = np.dot(inbuf, adjusted_weights) + adjusted_biases

    def _backwardImplementation(self, outerr, inerr, inbuf):
        adjusted_weights, _ = self._adjust_weights_biases(len(inbuf))
        inerr[:] = np.dot(outerr, adjusted_weights.T)

min_size = 5
max_size = 20
layer_size = 10

custom_layer = CustomNeuronLayer(layer_size, min_size, max_size)

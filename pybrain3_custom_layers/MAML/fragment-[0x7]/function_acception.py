# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )


"""

Accepts functions to be used within a list for each forward and backward pass. 

"""



from pybrain3.structure.modules.neuronlayer import NeuronLayer
import numpy as np

class CustomHiddenLayer(NeuronLayer):
		def __init__(self, size, functions=None, **kwargs):
				super(CustomHiddenLayer, self).__init__(size, **kwargs)
				self.functions = functions if functions is not None else []
				self.attention_weights = np.ones(size)

		def _forwardImplementation(self, inbuf, outbuf):
				for func in self.functions:
						outbuf[:] = func(inbuf)
						scores = np.dot(inbuf, self.attention_weights)
						attention_weights = np.exp(scores) / np.sum(np.exp(scores))
				outbuf[:] = np.dot(inbuf, attention_weights)
			
		def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
				for func in self.functions:
						outbuf[:] = func(inbuf)
						grad_attention_weights = np.dot(inbuf.T, outerr)
						self.attention_weights -= 0.001 * grad_attention_weights

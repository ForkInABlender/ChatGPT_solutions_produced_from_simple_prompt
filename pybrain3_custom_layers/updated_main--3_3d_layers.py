# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
This is similar to the gpt_script.py and functions the exact same way but with far more training needed.

This is tiktoken and pybrain3.

Tiktoken does the encoding, numpy parts do the work for BPE of the response, then it gets decoded. With no change of input, the output
 can be expected to be the same. As it also needs a feed back loop, it probably is a good idea to use recurrent or concurrent networks,
  or both. All depends on use-case.

Never-the-less this is another example of how gpt-2 functions at its core outside its normal housing.

Because it is using lightweight mechanical modeling for static computational graphs for something dynamic, it must also pay attention to
 to its own formulation of response.

From here, it is a matter of finding a dataset that doesn't use fake or stolen data/content.
"""

from pybrain3.structure import LinearLayer, RecurrentNetwork, FeedForwardNetwork, FullConnection
from dim3_neuronlayer import Dim3NeuronLayer  # Ensure this matches your file structure

import pickle
import tiktoken

class DynamicNet(RecurrentNetwork):
    def __init__(self, vocab_size, hidden_dim, num_heads=10):  # Added num_heads parameter with a default value
        super(DynamicNet, self).__init__()
        # Layers
        self.inLayer = LinearLayer(vocab_size)
        
        self.dim3Layer1 = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="d3_layer1")
        self.dim3Layer = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="d3_layer")  # Custom 3D layer
        self.dim3Layer3 = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="d3_layer3")
        self.outLayer = LinearLayer(vocab_size)
        
        # Adding modules and connections
        self.addInputModule(self.inLayer)
        
        self.addModule(self.dim3Layer1)
        self.addModule(self.dim3Layer)  # Add the custom 3D layer to the network
        self.addModule(self.dim3Layer3)
        
        self.addOutputModule(self.outLayer)
        
        self.in_to_dim3 = FullConnection(self.inLayer, self.dim3Layer1)  # Connection to the custom 3D layer
        
        self.dim3o1_to_dim3_h = FullConnection(self.dim3Layer1, self.dim3Layer)
        self.dim3_h_to_dim3o3 = FullConnection(self.dim3Layer, self.dim3Layer3)
        
        self.dim3_to_out = FullConnection(self.dim3Layer3, self.outLayer)
        
        self.addConnection(self.in_to_dim3)
        self.addConnection(self.dim3o1_to_dim3_h)
        self.addRecurrentConnection(FullConnection(self.dim3Layer, self.dim3Layer))
        self.addConnection(self.dim3_h_to_dim3o3)
        self.addConnection(self.dim3_to_out)
        
        self.sortModules()

def generate_response(model, input_text, tokenizer):
    input_tokens = tokenizer.encode(input_text)
    input_vector = [0] * vocab_size
    for token in input_tokens:
        input_vector[token] = 1
    output_tokens = model.activate(input_vector)
    return tokenizer.decode([int(output_tokens.argmax())])  # Ensure decoding is compatible with your tokenizer's API

if __name__ == "__main__":
    tokenizer = tiktoken.encoding_for_model("gpt-4")
    vocab_size = 50256
    hidden_dim = 128
    num_heads = 16 # Example number of heads for the custom 3D layer
    model = DynamicNet(vocab_size, hidden_dim, num_heads)
    input_text = "Hello, how are you?"
    response = generate_response(model, input_text, tokenizer)
    print("Response:", response)
    with open('dynamic_net.pkl', 'wb') as file:
        pickle.dump(model, file)
    print("Network configuration saved as JSON string to dynamic_net.pkl")
    print("saving complete")
    
    
    """
To unpickle the model:

``
import pickle

# Load the model from disk
with open('dynamic_net.pkl', 'rb') as file:
    loaded_model = pickle.load(file)

print("Model loaded from 'dynamic_net.pkl'")
``

    """

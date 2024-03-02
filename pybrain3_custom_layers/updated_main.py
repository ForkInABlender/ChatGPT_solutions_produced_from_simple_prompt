# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
Now it is updated with a neural network layer incorporating functions to emulate neurological responses chemically as well as their activation states.


"""

from pybrain3.structure import LinearLayer, RecurrentNetwork, FeedForwardNetwork, FullConnection
from dim3_neuronlayer import Dim3NeuronLayer  # Ensure this matches your file structure
from MolecularNeuroLayer import MolecularNeuroModule
from pybrain3.tools.xml.networkwriter import NetworkWriter

import tiktoken

class DynamicNet(RecurrentNetwork):
    def __init__(self, vocab_size, hidden_dim, num_heads=10):  # Added num_heads parameter with a default value
        super(DynamicNet, self).__init__()
        # Layers
        self.inLayer = LinearLayer(vocab_size)
        
        self.dim3Layer1 = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="d3_layer1")
        self.caffinespiderdance = MolecularNeuroModule(hidden_dim, hidden_dim, name="nice_web")
        self.dim3Layer = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="d3_layer")  # Custom 3D layer
        self.dim3Layer3 = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="d3_layer3")
        self.outLayer = LinearLayer(vocab_size)
        
        # Adding modules and connections
        self.addInputModule(self.inLayer)
        
        self.addModule(self.dim3Layer1)
        self.addModule(self.caffinespiderdance)
        self.addModule(self.dim3Layer)  # Add the custom 3D layer to the network
        self.addModule(self.dim3Layer3)
        
        self.addOutputModule(self.outLayer)
        
        self.in_to_dim3 = FullConnection(self.inLayer, self.dim3Layer1)  # Connection to the custom 3D layer
        self.dim3o1_to_niceweb = FullConnection(self.dim3Layer1, self.caffinespiderdance)
        self.dim3o1_to_dim3_h = FullConnection(self.caffinespiderdance, self.dim3Layer)
        self.dim3_h_to_dim3o3 = FullConnection(self.dim3Layer, self.dim3Layer3)
        
        self.dim3_to_out = FullConnection(self.dim3Layer3, self.outLayer)
        
        self.addConnection(self.in_to_dim3)
        self.addConnection(self.dim3o1_to_niceweb)
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
    vocab_size = 50257
    hidden_dim = 128
    num_heads = 16 # Example number of heads for the custom 3D layer
    model = DynamicNet(vocab_size, hidden_dim, num_heads)
    input_text = "Hello, how are you?"
    response = generate_response(model, input_text, tokenizer)
    print("Response:", response)
    NetworkWriter.writeToFile(model, '/app/dynamic_net.xml')
    print("Network configuration saved as xml file: '/app/dynamic_net.xml'")
    print("saving complete")

"""
To load the model from pybrain3's xml files:

``
from pybrain3.tools.customxml.networkreader import NetworkReader
model = NetworkReader.readFrom('/app/dynamic_net.xml')

``

"""

# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
This one exploded due too much complexity all at once being forced into the mold.

Instead I have to make 6 to 7 models that work as a larger model. Meaning I have to have the mirror neurons, prefrontal and frontal using a different gpt that uses
Dim3NeuronLayer, etcetera.



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
        # Define all custom layers with detailed configuration
        self.insula_net_layer = Dim3NeuronLayer(vocab_size, 120, num_heads=1, name="insula")
        self.frontal_lobe_net_layer = Dim3NeuronLayer(120, 100, num_heads=12, name="frontal_lobe")
        self.prefrontal_cortex_layer = Dim3NeuronLayer(250, 100, num_heads=10, name="prefrontal_lobe")
        self.cerebral_cortex_layer = Dim3NeuronLayer(925, 100, num_heads=5, name="cerebral_cortex")
        self.temporal_lobe_layer = Dim3NeuronLayer(50, 250, num_heads=25, name="temporal_lobe")
        self.visual_cortex_layer = Dim3NeuronLayer(120, 80, num_heads=8, name="visual_cortex")
        self.brocas_area_net_layer = Dim3NeuronLayer(150, 75, num_heads=5, name="brocas_area")
        self.wernickes_area_net_layer = Dim3NeuronLayer(75, 70, num_heads=5, name="wernickes_area")
        self.cerebellum_layer = Dim3NeuronLayer(700, 70, num_heads=7, name="cerebellum")
        self.brainstem_layer = Dim3NeuronLayer(120, 60, num_heads=6, name="brainstem")
        self.hippocampus_layer = Dim3NeuronLayer(250, 50, num_heads=5, name="hippocampus")
        self.thalamus_layer = Dim3NeuronLayer(400, 120, num_heads=10, name="thalamus")
        self.amygdala_layer = Dim3NeuronLayer(120, 10, num_heads=1, name="amygdala")
        self.hypothalamus_layer = Dim3NeuronLayer(120, 20, num_heads=2, name="hypothalamus")
        self.mirror_net_layer = Dim3NeuronLayer(500, 250, num_heads=25, name="mirror")
        self.outLayer = LinearLayer(vocab_size)
        # Adding modules and connections
        self.addInputModule(self.inLayer)
        self.addModule(self.insula_net_layer)
        self.addModule(self.frontal_lobe_net_layer)
        self.addModule(self.prefrontal_cortex_layer)
        self.addModule(self.cerebral_cortex_layer)
        self.addModule(self.temporal_lobe_layer)
        self.addModule(self.visual_cortex_layer)
        self.addModule(self.brocas_area_net_layer)
        self.addModule(self.wernickes_area_net_layer)
        self.addModule(self.cerebellum_layer)
        self.addModule(self.brainstem_layer)
        self.addModule(self.hippocampus_layer)
        self.addModule(self.thalamus_layer)
        self.addModule(self.amygdala_layer)
        self.addModule(self.hypothalamus_layer)
        self.addModule(self.mirror_net_layer)
        self.addOutputModule(self.outLayer)
        
        self.addConnection(FullConnection(self.inLayer, self.wernickes_area_net_layer))
        
        self.addConnection(FullConnection(self.inLayer, self.insula_net_layer))
        self.addConnection(FullConnection(self.insula_net_layer, self.frontal_lobe_net_layer))
        self.addConnection(FullConnection(self.frontal_lobe_net_layer, self.prefrontal_cortex_layer))
        self.addConnection(FullConnection(self.prefrontal_cortex_layer, self.cerebral_cortex_layer))
        self.addConnection(FullConnection(self.cerebral_cortex_layer, self.temporal_lobe_layer))
        self.addConnection(FullConnection(self.temporal_lobe_layer, self.visual_cortex_layer))
        self.addConnection(FullConnection(self.visual_cortex_layer, self.brocas_area_net_layer))
        self.addConnection(FullConnection(self.brocas_area_net_layer, self.wernickes_area_net_layer))
        self.addConnection(FullConnection(self.wernickes_area_net_layer, self.cerebellum_layer))
        self.addConnection(FullConnection(self.cerebellum_layer, self.brainstem_layer))
        self.addConnection(FullConnection(self.brainstem_layer, self.hippocampus_layer))
        self.addConnection(FullConnection(self.hippocampus_layer, self.thalamus_layer))
        self.addConnection(FullConnection(self.thalamus_layer, self.amygdala_layer))
        self.addConnection(FullConnection(self.amygdala_layer, self.hypothalamus_layer))
        self.addConnection(FullConnection(self.hypothalamus_layer, self.mirror_net_layer))
        self.addConnection(FullConnection(self.mirror_net_layer, self.outLayer))
        
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
    hidden_dim = 768
    num_heads = 16  # Example number of heads for the custom 3D layer
    model = DynamicNet(vocab_size, hidden_dim, num_heads)
    input_text = "Hello, how are you?"
    response = generate_response(model, input_text, tokenizer)
    print("Response:", response)
    with open('dynamic_net.pkl', 'wb') as file:
        pickle.dump(model, file)
    print("Network configuration saved as JSON string to dynamic_net.pkl")
    print("Saving complete")
    
    """
To unpickle the model:

import pickle

# Load the model from disk
with open('dynamic_net.pkl', 'rb') as file:
    loaded_model = pickle.load(file)

print("Model loaded from 'dynamic_net.pkl'")
    """

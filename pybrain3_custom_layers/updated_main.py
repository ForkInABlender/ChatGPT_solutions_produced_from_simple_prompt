# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
Now it is updated with a neural network layer incorporating functions to emulate neurological responses chemically as well as their activation states.

Now it can be trained on a large corpus of data, in practice, be useful for giving it, while emulated, the neurochemical equivalent to that of humans.

The last edits to make require edits to the ``MolecularNeuroModule`` class so one can assign them to assign specific chemicals used for signalling by the brain.

But it now also runs entirely offline as well.

"""


import os

tiktoken_cache_dir = "/app/"
os.environ["TIKTOKEN_CACHE_DIR"] = tiktoken_cache_dir
assert os.path.exists(os.path.join(tiktoken_cache_dir,"9b5ad71b2ce5302211f9c61530b329a4922fc6a4"))


from pybrain3.structure import LinearLayer, SigmoidLayer, TanhLayer, RecurrentNetwork, FeedForwardNetwork, FullConnection
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer
from dim3_neuronlayer import Dim3NeuronLayer
from numba import prange
import numpy as np
import tiktoken
import pickle

def pad_encoded_sequence(encoded_sequence, max_length, pad_value=0):
    num_padding = max_length - len(encoded_sequence)
    padded_sequence = encoded_sequence + [pad_value] * num_padding
    return np.array(padded_sequence, dtype=np.int32)

def generate_response(model, input_text, tokenizer):
    input_tokens = tokenizer.encode(input_text)
    input_vector = pad_encoded_sequence(input_tokens, vocab_size)
    output_tokens = model.activate(input_vector)
    return tokenizer.decode([int(output_tokens.argmax())]), [output_tokens]


class DynamicNet(RecurrentNetwork):
    def __init__(self, vocab_size, hidden_dim, num_heads=10):  # Added num_heads parameter with a default value
        super(DynamicNet, self).__init__()
        # Layers
        self.inLayer = LinearLayer(vocab_size)
        self.dim3Layer = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads, name="3d")  # Custom 3D layer
        self.outLayer = LinearLayer(vocab_size)
        # Adding modules and connections
        self.addInputModule(self.inLayer)
        self.addModule(self.dim3Layer)  # Add the custom 3D layer to the network
        self.addOutputModule(self.outLayer)
        self.in_to_dim3 = FullConnection(self.inLayer, self.dim3Layer)  # Connection to the custom 3D layer
        self.dim3_to_out = FullConnection(self.dim3Layer, self.outLayer)
        self.addConnection(self.in_to_dim3)
        self.addRecurrentConnection(FullConnection(self.dim3Layer, self.dim3Layer))
        self.addConnection(self.dim3_to_out)
        self.sortModules()

tokenizer = tiktoken.encoding_for_model("gpt-4")
vocab_size = 50257
hidden_dim = 29
num_heads = 1 # 1 head per hidden dimension squared. 
scale=vocab_size**3
model = DynamicNet(vocab_size, hidden_dim, num_heads)
model.randomize()
ds = SupervisedDataSet(vocab_size, vocab_size)
in_enc, out_enc = tokenizer.encode("Hello, how are you?"), tokenizer.encode("Fine. How are you?")
ds.addSample(pad_encoded_sequence(in_enc, vocab_size)/scale, pad_encoded_sequence(out_enc, vocab_size)/scale)

trainer = BackpropTrainer(model, ds)
for epoch in prange(1000000):
    print(epoch, trainer.train(), generate_response(model, "Hello, how are you?", tokenizer))
    
pickle.dump(model, open("/app/dynamic_net.pickle", "wb"))
print("saving complete")
print("done")


"""
To load the model from pickle files:

``
import pickle
model = pickle.load(open("/app/dynamic_net.pickle", "rb"))
``


"""

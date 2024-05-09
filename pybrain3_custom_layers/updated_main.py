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


from pybrain3.structure import LinearLayer, SigmoidLayer, TanhLayer, GaussianLayer, LSTMLayer, RecurrentNetwork, FeedForwardNetwork, FullConnection
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer
from main__GPT_Model import EmbeddingLayer
from dim3_neuronlayer import Dim3NeuronLayer
from numba import prange
import numpy as np
import tiktoken
import pickle

def pad_encoded_sequence(encoded_sequence, max_length, pad_value=220):
    num_padding = max_length - len(encoded_sequence)
    padded_sequence = encoded_sequence + [pad_value] * num_padding
    return np.array(padded_sequence, dtype=np.int64)

def generate_response(model, input_text, tokenizer):
    input_tokens = tokenizer.encode(input_text)
    input_vector = pad_encoded_sequence(input_tokens, vocab_size)
    output_tokens = model.activate(input_vector)
    return tokenizer.decode([abs(int(round(a.tolist())))//100000 for a in output_tokens]) # This change allows it to actually give more than one word responses....
#                                                          ^
#                                                          |
# Increase this number as the number of heads within the hidden dimension are used. For instance, 768 hidden dimensions in the 3-d layer will mean each node is treated as an
#  attention head while pertaining to 3-d conceptualization of information or otherwise to a value of 100000. 
# Modifying either the hidden_dim value and this number can also tell you how much memory & CPU is used. 
#
# If you are on linux or WSL and using htop, you'll be able to see how much it is using. On my host desktop machine, it runs slower due to this it-factor of how it runs.
#  Because of how small it is, overtime, it should be efficient enough to be useful after training has taken place. Meaning it can be run inside of a google cloud enviornment
#   or run right on GKE or AWS and the like using docker.
#
#
#
#

class DynamicNet(RecurrentNetwork):
    def __init__(self, vocab_size, hidden_dim, num_heads):
        super(DynamicNet, self).__init__()
        # Define the layers
        self.inLayer = LinearLayer(vocab_size)
        #self.embeddingLayer = embedding_layer  # Add the embedding layer
        self.dim3Layer = Dim3NeuronLayer(hidden_dim, hidden_dim, num_heads=num_heads, name="3d layer")
        self.outLayer = LinearLayer(vocab_size)
        # Adding modules and connections
        self.addInputModule(self.inLayer)
        self.addModule(self.embeddingLayer)  # Add the embedding layer to the network
        self.addModule(self.dim3Layer)
        self.addOutputModule(self.outLayer)
        self.in_to_embedding = FullConnection(self.inLayer, self.dim3Layer)  # Connection to the embedding layer
        #self.embedding_to_dim3 = FullConnection(self.embeddingLayer, self.dim3Layer)  # Connection from embedding to custom 3D layer
        self.addConnection(self.in_to_embedding)
        #self.addConnection(self.embedding_to_dim3)
        for a in range(num_heads):
            self.addRecurrentConnection(FullConnection(self.dim3Layer, self.dim3Layer))
        self.dim3_to_out = FullConnection(self.dim3Layer, self.outLayer)
        self.addConnection(self.dim3_to_out)
        self.sortModules()

tokenizer = tiktoken.encoding_for_model("gpt-4")
vocab_size = 100257
hidden_dim = 768
num_heads = 16 # number of heads per hidden dimension cubed. Or roughly 452 billion neurons for 768 heads (all heads) allotted
scale=int(vocab_size**3.61)
ds = SupervisedDataSet(vocab_size, vocab_size)
#embedding_layer = embedding_layer = EmbeddingLayer(vocab_size, 128)
model = DynamicNet(vocab_size, hidden_dim, num_heads) ## model = pickle.load(open("/app/dynamic_net.pickle", "rb")) # for loading from a pickle file.
trainer = BackpropTrainer(model, ds, learningrate=.0931) # Give it a 9.31% learning rate to keep it from committing to over-fitting.

def gen_qa(user, machine): # develop Q&A for the response type it should have. Including function calls. This means that function calls will need to watch for 'function_call',
    in_enc, out_enc = tokenizer.encode(user), tokenizer.encode(machine) # it's name, and parameters.
    ds.addSample(pad_encoded_sequence(in_enc, vocab_size)/scale, pad_encoded_sequence(out_enc, vocab_size)/scale)

for epoch in range(150000): # Give it time to formulate response after training by the epoch. This lets one gage where the curve for it is without jumping the gun.
    print(f"epoch: {epoch}\ntrain: {trainer.train()}\nResponse: {generate_response(model, 'Hello, how are you?', tokenizer)}")

pickle.dump(model, open("/app/dynamic_net.pickle", "wb")) # to save the model after training. 
print("saving complete")
print("done")

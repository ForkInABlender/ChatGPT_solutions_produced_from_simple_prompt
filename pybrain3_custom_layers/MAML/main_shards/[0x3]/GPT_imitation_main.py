# Dylan Kenneth Eliot

"""
This creates with the Dim3NeuronLayer as the hidden layers. This should form the basic template for GPT knowing that the Dim3NeuronLayer does the QKV as well as 3d, 2d, & 1d processing of data.

The Dim3NeuronLayer itself acts as a layer for allowing a 3 dimensional neuron comprised of its 2 dimensional counter parts, using them as anchor corner points.
 This gives the FFT the ability to update the connections and individual data points, allowing for 3d processing, spatial awareness, image processing, and language processing.

language datasets soon to follow as well as brain-emulator code revision.


update 10/19/2024 -- @ 3:01:03 pm -- This file is now for use in AI modeling purposes. For now, this is for imitating GPT-3.5, 4, 4o, & 4+. Or at the very least, the parts that are core to main brain function of the stoma.

update 10/22/2024 -- 2:51:21 pm -- All parts now function as expected.

"""


import numpy as np
from pybrain3.structure import FeedForwardNetwork, LinearLayer, FullConnection, Module
from pybrain3.structure.modules.neuronlayer import NeuronLayer
from dim3_neuronlayer import Dim3NeuronLayer

def softmax(x, axis=-1):
    e_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
    return e_x / e_x.sum(axis=axis, keepdims=True)

class CustomBatchLayer(Module):
    """Custom layer to handle batch input with _forwardImplementation and _backwardImplementation."""

    def __init__(self, dim_in, dim_out):
        super(CustomBatchLayer, self).__init__(dim_in, dim_out)

    def _forwardImplementation(self, inbuf, outbuf):
        batch_size = inbuf.size // self.indim
        input_matrix = np.reshape(inbuf, (batch_size, self.indim))
        activated_matrix = 1 / (1 + np.exp(-np.clip(input_matrix, -700, 700)))  # Sigmoid with clipping
        np.copyto(outbuf, activated_matrix.flatten())

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        batch_size = inbuf.size // self.indim
        outbuf_matrix = np.reshape(outbuf, (batch_size, self.indim))
        sigmoid_grad = outbuf_matrix * (1 - outbuf_matrix)
        np.copyto(inerr, outerr * sigmoid_grad.flatten())

class EmbeddingLayer(LinearLayer):
    def __init__(self, vocab_size, embedding_dim):
        super().__init__(embedding_dim)
        self.embeddings = np.random.randn(vocab_size, embedding_dim)

    def _forwardImplementation(self, inbuf, outbuf):
        self.token_idx = np.argmax(inbuf)
        outbuf[:] = self.embeddings[min(self.token_idx, len(self.embeddings) - 1)]

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        gradient = np.zeros_like(self.embeddings)
        gradient[self.token_idx] = outerr
        self.embeddings -= gradient
        inerr[:] = self.embeddings.T @ outerr


class MultiHeadSelfAttention(NeuronLayer):
    def __init__(self, indim, outdim, num_heads):
        super().__init__(indim, outdim)
        self.num_heads = num_heads
        self.depth = indim // num_heads
        self.W_q = np.random.rand(indim, indim)
        self.W_k = np.random.rand(indim, indim)
        self.W_v = np.random.rand(indim, indim)
        self.W_o = np.random.rand(indim, outdim)

    def scaled_dot_product_attention(self, Q, K, V):
        matmul_qk = np.dot(Q, K.T)
        d_k = Q.shape[-1]
        scaled_attention_logits = matmul_qk / np.sqrt(d_k)
        attention_weights = softmax(scaled_attention_logits, axis=-1)
        return np.dot(attention_weights, V)

    def _forwardImplementation(self, inbuf, outbuf):
        if len(inbuf.shape) == 1:
            inbuf = inbuf[np.newaxis, :]
        Q = np.dot(inbuf, self.W_q)
        K = np.dot(inbuf, self.W_k)
        V = np.dot(inbuf, self.W_v)
        Q = np.split(Q, self.num_heads, axis=1)
        K = np.split(K, self.num_heads, axis=1)
        V = np.split(V, self.num_heads, axis=1)

        attention_heads = []
        for i in range(self.num_heads):
            attention_head = self.scaled_dot_product_attention(Q[i], K[i], V[i])
            attention_heads.append(attention_head)

        # Concatenate attention heads and pass through the final linear layer
        concatenated_heads = np.concatenate(attention_heads, axis=1)
        outbuf[:] = np.dot(concatenated_heads, self.W_o)

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        d_concat_heads = np.dot(outerr, self.W_o.T)
        d_attention_heads = np.split(d_concat_heads, self.num_heads, axis=1)
        dQ_total, dK_total, dV_total = 0, 0, 0
        for i in range(self.num_heads):
            d_out = d_attention_heads[i]
            dV = np.dot(self.attention_weights[i].T, d_out)
            d_attention_weights = np.dot(d_out, self.V[i].T)
            dQK = d_attention_weights * (1 - self.attention_weights[i]) * self.attention_weights[i]
            dQ = np.dot(dQK, self.K[i])
            dK = np.dot(dQK, self.Q[i])
            dQ_total += dQ
            dK_total += dK
            dV_total += dV
        inerr[:] = np.dot(dQ_total, self.W_q.T) + np.dot(dK_total, self.W_k.T) + np.dot(dV_total, self.W_v.T)
        self.W_q -= np.dot(inbuf.T, dQ_total)
        self.W_k -= np.dot(inbuf.T, dK_total)
        self.W_v -= np.dot(inbuf.T, dV_total)
        self.W_o -= np.dot(self.concatenated_heads.T, outerr)

class LayerNorm(NeuronLayer):
    def __init__(self, dim, eps=1e-5):
        super().__init__(dim, dim)
        self.eps = eps
        self.gamma = np.ones(dim)
        self.beta = np.zeros(dim)

    def _forwardImplementation(self, inbuf, outbuf):
        mean = np.mean(inbuf, axis=-1, keepdims=True)
        variance = np.var(inbuf, axis=-1, keepdims=True)
        outbuf[:] = self.gamma * (inbuf - mean) / np.sqrt(variance + self.eps) + self.beta

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        # Placeholder backward implementation
        inerr[:] = outerr  # In real cases, more sophisticated differentiation is needed.


# Constants
VOCAB_SIZE = 50256
D_MODEL = 128
FFN_DIM = 128
NUM_HEADS = 128
NUM_BLOCKS = 2100

# Create the GPT-like network
net = FeedForwardNetwork()
inLayer = LinearLayer(VOCAB_SIZE)
net.addInputModule(inLayer)

# Add embedding layer and custom batch layer
embedding = EmbeddingLayer(VOCAB_SIZE, D_MODEL)
batch_layer = CustomBatchLayer(D_MODEL, D_MODEL)
net.addModule(embedding)
net.addModule(batch_layer)
net.addConnection(FullConnection(inLayer, embedding))
net.addConnection(FullConnection(embedding, batch_layer))

# Add multi-head attention and transformer blocks
prev_layer = batch_layer
for _ in range(NUM_BLOCKS):
    norm1 = LayerNorm(D_MODEL)
    attention = MultiHeadSelfAttention(D_MODEL, D_MODEL, NUM_HEADS)
    dim3_layer = Dim3NeuronLayer(
        in_dim=D_MODEL, out_dim=D_MODEL, embed_dropout=0.1, attn_dropout=0.1,
        activation_function='gelu', lr=0.0001, weight_decay=0.01, gradient_clipping=1.0
    )
    norm2 = LayerNorm(D_MODEL)
    
    net.addModule(norm1)
    net.addModule(attention)
    net.addModule(dim3_layer)
    net.addModule(norm2)

    net.addConnection(FullConnection(prev_layer, norm1))
    net.addConnection(FullConnection(norm1, attention))
    net.addConnection(FullConnection(attention, dim3_layer))
    net.addConnection(FullConnection(dim3_layer, norm2))
    
    prev_layer = norm2

# Add final output layer
out_layer = LinearLayer(VOCAB_SIZE)
net.addOutputModule(out_layer)
net.addConnection(FullConnection(prev_layer, out_layer))

# Finalize the network
net.sortModules()





# Pad to cap-length (VOCAB_Size); treat `-1` as filler for word delimination on 'encoding' from ordinal sets of each word, splitting on space only.
def pad_to_length(array, target_length=VOCAB_SIZE, pad_value=-1):
    """Pad each word's ordinal array to the target length."""
    return array + [pad_value] * (target_length - len(array))


# Example padded ordinal array
padded_ordinals = [] # put up to the vocab limit, split by space character, for proper padding before use.
# Flatten the batch input (required for network activation)
flat_input = np.array(pad_to_length(padded_ordinals))

# Activate the network with the flattened input
output = net.activate(flat_input)

# Print the network's output
print("Network output for the input batch:")
print(output)

#print(reshaped_output)

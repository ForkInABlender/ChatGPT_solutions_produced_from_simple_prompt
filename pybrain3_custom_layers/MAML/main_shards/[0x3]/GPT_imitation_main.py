# Dylan Kenneth Eliot

"""
This creates with the Dim3NeuronLayer as the hidden layers. This should form the basic template for GPT knowing that the Dim3NeuronLayer does the QKV as well as 3d, 2d, & 1d processing of data.

The Dim3NeuronLayer itself acts as a layer for allowing a 3 dimensional neuron comprised of its 2 dimensional counter parts, using them as anchor corner points.
 This gives the FFT the ability to update the connections and individual data points, allowing for 3d processing, spatial awareness, image processing, and language processing.

language datasets soon to follow as well as brain-emulator code revision.


update 10/19/2024 -- @ 3:01:03 pm -- This file is now for use in AI modeling purposes. For now, this is for imitating GPT-3.5, 4, 4o, & 4+. Or at the very least, the parts that are core to main brain function of the stoma.
"""

from pybrain3.structure import FeedForwardNetwork, LinearLayer, FullConnection
from pybrain3.structure import SoftmaxLayer
from pybrain3.structure.modules.neuronlayer import NeuronLayer
from dim3_neuronlayer import Dim3NeuronLayer
import numpy as np

class EmbeddingLayer(LinearLayer):
    def __init__(self, vocab_size, embedding_dim):
        super().__init__(embedding_dim)
        self.embeddings = np.random.randn(vocab_size, embedding_dim)

    def _forwardImplementation(self, inbuf, outbuf):
        self.token_idx = np.argmax(inbuf)
        if self.token_idx <= len(self.embeddings):
            outbuf[:] = self.embeddings[self.token_idx]
        else:
            self.token_idx = len(self.embeddings) - 1
            outbuf[:] = self.embeddings[self.token_idx]

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        gradient = np.zeros_like(self.embeddings)
        gradient[self.token_idx] = outerr
        self.embeddings -= gradient
        inerr[:] = self.embeddings.T @ outerr

class GeLULayer(NeuronLayer):
    def __init__(self, dim):
        super().__init__(dim, dim)

    def _forwardImplementation(self, inbuf, outbuf):
        outbuf[:] = 0.5 * inbuf * (
            1 + np.tanh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf ** 3)))

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        cdf = 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf ** 3)))
        pd = np.sqrt(2 / np.pi) * (1 + 0.134145 * inbuf ** 2) * (1 / np.cosh(np.sqrt(2 / np.pi) * (inbuf + 0.044715 * inbuf ** 3))) ** 2
        inerr[:] = outerr * (cdf + inbuf * pd)

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

# Create the GPT-like network with updated architecture
VOCAB_SIZE = 50257  # Typical GPT vocabulary size
D_MODEL = 128  # Embedding size for GPT-like model; median for x86_64 is 768
FFN_DIM = 128  # Feed-forward network dimension in GPT blocks
NUM_HEADS = 12
NUM_BLOCKS = 96

net = FeedForwardNetwork()
inLayer = LinearLayer(VOCAB_SIZE)
net.addInputModule(inLayer)

# Embedding layer to convert input tokens to vector representations
embedding = EmbeddingLayer(VOCAB_SIZE, D_MODEL)
net.addModule(embedding)
net.addConnection(FullConnection(inLayer, embedding))

# Multi-head attention layer
attention = MultiHeadSelfAttention(D_MODEL, D_MODEL, NUM_HEADS)
net.addModule(attention)
net.addConnection(FullConnection(embedding, attention))

prev_layer = attention

# Add transformer blocks with Dim3NeuronLayer
for _ in range(NUM_BLOCKS):
    # Layer Normalization before Attention
    norm1 = LayerNorm(D_MODEL)
    net.addModule(norm1)
    net.addConnection(FullConnection(prev_layer, norm1))

    # Dim3NeuronLayer for multi-dimensional input
    dim3_layer = Dim3NeuronLayer(in_dim=D_MODEL, out_dim=D_MODEL, embed_dropout=0.1, attn_dropout=0.1, activation_function='gelu', lr=0.0001, weight_decay=0.01, gradient_clipping=1.0)
    net.addModule(dim3_layer)
    net.addConnection(FullConnection(norm1, dim3_layer))

    # Feed-Forward Layer 1
    ffn1 = LinearLayer(FFN_DIM, FFN_DIM)
    net.addModule(ffn1)
    net.addConnection(FullConnection(dim3_layer, ffn1))

    # GeLU Activation
    gelu = GeLULayer(FFN_DIM)
    net.addModule(gelu)
    net.addConnection(FullConnection(ffn1, gelu))

    # Feed-Forward Layer 2
    ffn2 = LinearLayer(FFN_DIM, FFN_DIM)
    net.addModule(ffn2)
    net.addConnection(FullConnection(gelu, ffn2))

    # Layer Normalization after Feed-Forward Network
    norm2 = LayerNorm(FFN_DIM)
    net.addModule(norm2)
    net.addConnection(FullConnection(ffn2, norm2))
    
    # Update previous layer to point to the latest norm layer
    prev_layer = norm2

# Add final output layer
out_layer = LinearLayer(D_MODEL)
net.addOutputModule(out_layer)
final_connection = FullConnection(prev_layer, out_layer)
net.addConnection(final_connection)

net.sortModules()

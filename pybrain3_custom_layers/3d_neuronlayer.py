# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""

What this allows for is neural networks to make use of 3-d information.

This is useful for when one is modeling the hippocampus, temporal, or paradial lobes as well as the visual cortex.

This is also meant for allowing lightweight LLMs using pybrain3 to make use of visual, auditorial, or sensory data essential to having an imagination or free thinking.

This also means that it must have a 3-dimensional modeling of even the cerebral cortex.


Sadly, this isn't meant to replace modern-day GPUs. However, if you find a use for it in that domain, more power to you. But it is to make use of the CPU as it "thinks" on a 
 thing it is looking at, listening to, or trying to reason about, when it comes to AI development. (Alan Turings' definition is the only acceptable definition. tough, cope. 

Token parsing is not AI.  And any failure to accept that will not be treated with respect as the Eliza projects on the market aren't AI. )

For more information, please read through

https://chat.openai.com/share/837e66b0-6ff7-458b-9489-5eb950449ae8

Enjoy :)

You will not be missed!

"""

import numpy as np
from pybrain.structure.modules import NeuronLayer

class Dim3NeuronLayer(NeuronLayer):
    def __init__(self, indim, outdim, num_heads):
        assert indim % num_heads == 0, "indim must be divisible by num_heads"
        super(Custom3DNeuronLayer, self).__init__()
        self.indim = indim
        self.outdim = outdim
        self.num_heads = num_heads
        self.depth = indim // num_heads

        self.W_q = np.random.randn(indim, indim)
        self.W_k = np.random.randn(indim, indim)
        self.W_v = np.random.randn(indim, indim)
        self.W_o = np.random.randn(indim, outdim)

    def scaled_dot_product_attention(self, Q, K, V):
        matmul_qk = np.dot(Q, K.swapaxes(-2, -1))
        dk = Q.shape[-1]
        scaled_attention_logits = matmul_qk / np.sqrt(dk)
        
        # Softmax is applied to the last axis (seq_length_k) to normalize the scores
        attention_weights = np.exp(scaled_attention_logits) / np.sum(np.exp(scaled_attention_logits), axis=-1, keepdims=True)
        output = np.dot(attention_weights, V)
        return output, attention_weights

    def _forwardImplementation(self, inbuf, outbuf):
        batch_size, seq_length, _ = inbuf.shape

        # Linear projections in batch from input_dim => num_heads x depth
        Q = np.dot(inbuf.reshape(-1, self.indim), self.W_q).reshape(batch_size, seq_length, self.num_heads, self.depth)
        K = np.dot(inbuf.reshape(-1, self.indim), self.W_k).reshape(batch_size, seq_length, self.num_heads, self.depth)
        V = np.dot(inbuf.reshape(-1, self.indim), self.W_v).reshape(batch_size, seq_length, self.num_heads, self.depth)

        # Transpose to get dimensions batch_size x num_heads x seq_length x depth
        Q = Q.transpose(0, 2, 1, 3)
        K = K.transpose(0, 2, 1, 3)
        V = V.transpose(0, 2, 1, 3)

        # Scaled Dot-Product Attention
        attention_outputs, self.attention_weights = [], []
        for i in range(self.num_heads):
            attention_output, attention_weight = self.scaled_dot_product_attention(Q[:, i, :, :], K[:, i, :, :], V[:, i, :, :])
            attention_outputs.append(attention_output)
            self.attention_weights.append(attention_weight)

        # Concatenation of heads
        attention_output = np.concatenate(attention_outputs, axis=-1)

        # Final linear layer
        outbuf[:] = np.dot(attention_output.reshape(batch_size, -1), self.W_o)

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        batch_size, seq_length, _ = inbuf.shape

        # Gradients for final linear layer
        dW_o = np.dot(outbuf.reshape(-1, self.outdim).T, outerr)
        doutbuf = np.dot(outerr, self.W_o.T).reshape(batch_size, seq_length, self.num_heads, self.depth)

        # Split doutbuf for heads
        doutbuf = doutbuf.transpose(0, 2, 1, 3)  # Transpose back to num_heads first for consistency

        dW_q, dW_k, dW_v = np.zeros_like(self.W_q), np.zeros_like(self.W_k), np.zeros_like(self.W_v)
        for i in range(self.num_heads):
            dattention_output = doutbuf[:, i, :, :]

            # Placeholder for the gradient through the softmax and dot products
            # This is complex and involves partial derivatives through the attention mechanism
            # Assuming we have dQ, dK, dV from the attention gradients for simplification

            # Update gradients for Q, K, V weights
            dW_q_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))
            dW_k_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))
            dW_v_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))

            dW_q += dW_q_partial
            dW_k += dW_k_partial
            dW_v += dW_v_partial

        # Update weights (simplified without learning rate or optimization algorithm)
        self.W_q -= dW_q
        self.W_k -= dW_k
        self.W_v -= dW_v
        self.W_o -= dW_o

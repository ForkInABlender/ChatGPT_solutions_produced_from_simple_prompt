# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition) & Google-Gemini

"""

Now that the 3-d Layer works properly,

One can begin 3-d compute for each region with that needs it.

This was fixed due to an error I failed to spot. During training, this network now works as well.


* update note: had to be fixed for training and saving purposes. Reversion of optimizations for now. Apologies to those relying on optimizations. This is as optimized as it
 is likely to get. 

10/15/2024 -- This file is here for demonstration with the prior class within the same subdirectory { [0x3] }. Nothing has changed about this file.
"""




import numpy as np
from pybrain3.structure.modules import Module

class Dim3NeuronLayer(Module):
    def __init__(self, indim, outdim, num_heads, name=None):
        super(Dim3NeuronLayer, self).__init__(indim, outdim, name=name)
        self.indim = indim
        self.outdim = outdim
        self.num_heads = num_heads
        assert indim % num_heads == 0, "Input dimension must be divisible by the number of heads"
        self.depth = indim // num_heads
        self.seq_length = 1

        # Initialize weights for query, key, and value transformations
        self.W_q = np.random.randn(self.indim, self.indim)
        self.W_k = np.random.randn(self.indim, self.indim)
        self.W_v = np.random.randn(self.indim, self.indim)
        self.W_o = np.random.randn(self.depth * self.num_heads * self.seq_length, self.outdim)

        # Set input and output sizes
        self.inputbuffer = np.zeros((1, indim * self.seq_length))
        self.outputbuffer = np.zeros((1, outdim))

    def _forwardImplementation(self, inbuf, outbuf):
        if inbuf.ndim == 1:
            inbuf = inbuf.reshape(1, 1, -1)
        elif inbuf.ndim == 2:
            inbuf = inbuf.reshape(inbuf.shape[0], 1, -1)
        batch_size = inbuf.shape[0]

        # Reshape input to [batch_size, seq_length, indim]
        if inbuf.size != batch_size * self.seq_length * self.indim:
            raise ValueError(f"Input size mismatch: expected {batch_size * self.seq_length * self.indim}, got {inbuf.size}")

        inbuf_reshaped = inbuf.reshape(batch_size, self.seq_length, self.indim)

        # Compute Q, K, V matrices
        Q = np.dot(inbuf_reshaped, self.W_q).reshape(batch_size, self.seq_length, self.num_heads, self.depth)
        K = np.dot(inbuf_reshaped, self.W_k).reshape(batch_size, self.seq_length, self.num_heads, self.depth)
        V = np.dot(inbuf_reshaped, self.W_v).reshape(batch_size, self.seq_length, self.num_heads, self.depth)

        # Transpose to [batch_size, num_heads, seq_length, depth]
        self.Q = Q.transpose(0, 2, 1, 3)
        self.K = K.transpose(0, 2, 1, 3)
        self.V = V.transpose(0, 2, 1, 3)

        # Scaled dot-product attention for each head
        attention_outputs = []
        for i in range(self.num_heads):
            dk = self.depth
            matmul_qk = np.dot(self.Q[:, i, :, :], self.K[:, i, :, :].transpose(0, 2, 1))
            scaled_attention_logits = matmul_qk / np.sqrt(dk)
            attention_weights = np.exp(scaled_attention_logits - np.max(scaled_attention_logits, axis=-1, keepdims=True))
            attention_weights /= np.sum(attention_weights, axis=-1, keepdims=True)
            output = np.dot(attention_weights, self.V[:, i, :, :])
            attention_outputs.append(output)

        # Concatenate attention outputs from all heads
        self.attention_output = np.concatenate(attention_outputs, axis=-1)
        self.attention_output = self.attention_output.reshape(batch_size, -1)

        # Output projection
        outbuf[:] = np.dot(self.attention_output, self.W_o)


    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        batch_size = inbuf.shape[0] // self.seq_length

        # Reshape outerr to match the shape of attention_output
        outerr_reshaped = outerr.reshape(batch_size, self.depth * self.num_heads, self.seq_length)

        # Initialize gradients for weights
        dW_o = np.zeros_like(self.W_o)
        dW_q = np.zeros_like(self.W_q)
        dW_k = np.zeros_like(self.W_k)
        dW_v = np.zeros_like(self.W_v)

        # Gradient for output weights
        dW_o[:] = np.dot(self.attention_output.T, outerr)

        # Calculate gradients for the attention mechanism
        d_attention_output = np.dot(outerr, self.W_o.T).reshape(batch_size, self.num_heads, self.seq_length, self.depth)

        for i in range(self.num_heads):
            # Calculate gradients w.r.t. V
            d_attention_weights = np.dot(d_attention_output[:, i, :, :], self.V[:, i, :, :].transpose(0, 2, 1))
            dV = np.dot(self.Q[:, i, :, :].transpose(0, 2, 1), d_attention_weights)

            # Calculate gradients w.r.t. Q and K
            dQ = np.dot(d_attention_weights, self.K[:, i, :, :])
            dK = np.dot(d_attention_weights.transpose(0, 2, 1), self.Q[:, i, :, :])

            # Sum the gradients across batches
            dW_q += np.dot(inbuf.transpose(1, 0), dQ.reshape(batch_size * self.seq_length, self.depth))
            dW_k += np.dot(inbuf.transpose(1, 0), dK.reshape(batch_size * self.seq_length, self.depth))
            dW_v += np.dot(inbuf.transpose(1, 0), dV.reshape(batch_size * self.seq_length, self.depth))

        # Assign calculated gradients to inerr
        inerr[:] = np.dot(outerr, self.W_o.T)

        # Update weights with gradients (simple gradient descent)
        learning_rate = 0.001  # Example learning rate
        self.W_o -= learning_rate * dW_o
        self.W_q -= learning_rate * dW_q
        self.W_k -= learning_rate * dW_k
        self.W_v -= learning_rate * dW_v

# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition) & Google-Gemini

"""

Now that the 3-d Layer works properly,

One can begin 3-d compute for each region with that needs it.

This was fixed due to an error I failed to spot. During training, this network now works as well.


* update note: had to be fixed for training and saving purposes. Reversion of optimizations for now. Apologies to those relying on optimizations. This is as optimized as it
 is likely to get. 


"""

import numpy as np
from pybrain3.structure.modules.neuronlayer import NeuronLayer # forgot `.neuronlayer`.... oops
from concurrent.futures import ThreadPoolExecutor
from numba import prange

def softmax(x, axis=-1):
    e_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
    return e_x / e_x.sum(axis=axis, keepdims=True)

def compute_attention_for_head(args):
    Q, K, V, dk = args
    dk = max(dk, 1e-9)
    matmul_qk = np.dot(Q, K.swapaxes(-2, -1))
    matmul_qk = np.clip(matmul_qk, a_min=None, a_max=1e9)
    scaled_attention_logits = matmul_qk / np.sqrt(dk)
    logits_max = np.max(scaled_attention_logits, axis=-1, keepdims=True)
    logits_max = np.nan_to_num(logits_max, nan=0.0, posinf=1e9, neginf=-1e9)
    logits_diff = scaled_attention_logits - logits_max
    exp_logits = np.exp(np.clip(logits_diff, a_min=None, a_max=700))
    epsilon = 1e-9
    attention_weights = exp_logits / (np.sum(exp_logits, axis=-1, keepdims=True) + epsilon)
    output = np.dot(attention_weights, V)
    return output, attention_weights

class Dim3NeuronLayer(NeuronLayer):
    def __init__(self, indim, outdim, num_heads, name):
        assert indim % num_heads == 0, "indim must be divisible by num_heads"
        super(Dim3NeuronLayer, self).__init__(indim, outdim)
        self.indim = indim
        self.outdim = outdim
        self.num_heads = num_heads
        self.depth = indim // num_heads
        self.name=name
        self.W_q, self.W_k = np.random.randn(indim, indim)
        self.W_k = np.random.randn(indim, indim)
        self.W_v = np.random.randn(indim, indim)
        # Adjust the output dimension for concatenated heads
        self.W_o = np.random.randn(self.depth * num_heads, outdim)
    def scaled_dot_product_attention(self, Q, K, V):
        dk = Q.shape[-1]
        args_list = [(Q[:, i, :, :], K[:, i, :, :], V[:, i, :, :], dk) for i in range(Q.shape[1])] # Prepare arguments for parallel processing
        with ThreadPoolExecutor(max_workers=self.num_heads) as executor:  # Use multiprocessing pool to compute attention in parallel across heads # Adjust number of processes based on your environment
            results = list(executor.map(compute_attention_for_head, args_list))
        outputs, attention_weights = zip(*results)  # Separate outputs and attention weights from the results
        output = np.stack(outputs, axis=1) # Combine outputs from all heads
        return output, attention_weights  # Return both combined output and attention weights

    def _forwardImplementation(self, inbuf, outbuf):
        if inbuf.ndim == 1:
            inbuf = inbuf.reshape(1, 1, -1)
        elif inbuf.ndim == 2:
            inbuf = inbuf.reshape(inbuf.shape[0], 1, -1)
        batch_size, seq_length, _ = inbuf.shape
        Q = np.dot(inbuf.reshape(-1, self.indim), self.W_q).reshape(batch_size, seq_length, self.num_heads, self.depth)
        K = np.dot(inbuf.reshape(-1, self.indim), self.W_k).reshape(batch_size, seq_length, self.num_heads, self.depth)
        V = np.dot(inbuf.reshape(-1, self.indim), self.W_v).reshape(batch_size, seq_length, self.num_heads, self.depth)
        self.Q = Q.transpose(0, 2, 1, 3)
        self.K = K.transpose(0, 2, 1, 3)
        self.V = V.transpose(0, 2, 1, 3)
        attention_outputs, self.attention_weights = self.scaled_dot_product_attention(self.Q, self.K, self.V)
        self.attention_output = np.concatenate(attention_outputs, axis=-1)
        outbuf[:] = np.dot(self.attention_output.reshape(batch_size, -1), self.W_o)

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        if inbuf.ndim == 1:
            inbuf = inbuf.reshape(1, 1, -1)
        elif inbuf.ndim == 2:
            inbuf = inbuf.reshape(inbuf.shape[0], 1, -1)
        batch_size, seq_length, _ = inbuf.shape
        if outerr.ndim == 1:
            outerr = outerr.reshape(1, 1, -1)
        elif outerr.ndim == 2:
            outerr = outerr.reshape(outerr.shape[0], 1, -1)
        dW_o = np.dot(outbuf.reshape(-1, self.outdim).T, outerr.reshape(-1, self.outdim))
        dout_attention = np.dot(outerr, self.W_o.T).reshape(batch_size, self.num_heads, seq_length, self.depth)
        dout_attention = dout_attention.transpose(0, 2, 1, 3)
        dW_q = np.zeros_like(self.W_q)
        dW_k = np.zeros_like(self.W_k)
        dW_v = np.zeros_like(self.W_v)
        ##
        d_concat_heads = np.dot(outerr, self.W_o.T)
        d_attention_heads = np.split(d_concat_heads, self.num_heads, axis=2)
        dQ_total, dK_total, dV_total = 0, 0, 0
        for i in prange(self.num_heads):
            d_out = d_attention_heads[i]
            dV = np.dot(self.attention_weights[i].T, d_out)
            dV_total += dV
            d_attention_weights = np.dot(self.attention_weights, i)
            dQK = d_attention_weights * (1 - self.attention_weights[i]) * self.attention_weights[i]
            dQ = np.dot(dQK, self.K[:, i, :, :])
            dK = np.dot(dQK, self.Q[:, i, :, :])
            dV_total += dV
            dW_q += dQ.reshape(dW_q.shape[0])
            dW_k += dK.reshape(dW_k.shape[0])
        ##
        self.W_q -= dW_q
        self.W_k -= dW_k
        self.W_v -= dW_v
        self.W_o -= dW_o
        outbuf[:] = np.dot(self.attention_output.reshape(batch_size, -1), self.W_o)

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

def compute_attention_for_head(args):
    Q, K, V, dk = args
    matmul_qk = np.dot(Q, K.swapaxes(-2, -1))
    scaled_attention_logits = matmul_qk / np.sqrt(dk)
    logits_max = np.max(scaled_attention_logits, axis=-1, keepdims=True)
    exp_logits = np.exp(scaled_attention_logits - logits_max)
    attention_weights = exp_logits / np.sum(exp_logits, axis=-1, keepdims=True)
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

        # Initialize weights for the projections
        self.W_q = np.random.randn(indim, indim)
        self.W_k = np.random.randn(indim, indim)
        self.W_v = np.random.randn(indim, indim)
        # Adjust the output dimension for concatenated heads
        self.W_o = np.random.randn(self.depth * num_heads, outdim)

    def scaled_dot_product_attention(self, Q, K, V):
        dk = Q.shape[-1]
        # Prepare arguments for parallel processing
        args_list = [(Q[:, i, :, :], K[:, i, :, :], V[:, i, :, :], dk) for i in range(Q.shape[1])]
        
        # Use multiprocessing pool to compute attention in parallel across heads
        with ThreadPoolExecutor(max_workers=self.num_heads) as executor:  # Adjust number of processes based on your environment
            results = list(executor.map(compute_attention_for_head, args_list))
        
        # Separate outputs and attention weights from the results
        outputs, attention_weights = zip(*results)
        
        # Combine outputs from all heads
        output = np.stack(outputs, axis=1)
        
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

        Q = Q.transpose(0, 2, 1, 3)
        K = K.transpose(0, 2, 1, 3)
        V = V.transpose(0, 2, 1, 3)

        attention_outputs, attention_weights = self.scaled_dot_product_attention(Q, K, V)

        attention_output = np.concatenate(attention_outputs, axis=-1)
        outbuf[:] = np.dot(attention_output.reshape(batch_size, -1), self.W_o)

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

        for i in range(self.num_heads):
            dattention_output = dout_attention[:, :, i, :]
            dW_q_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))
            dW_k_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))
            dW_v_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))

            dW_q += -dW_q_partial
            dW_k += -dW_k_partial
            dW_v += -dW_v_partial

        self.W_q -= dW_q
        self.W_k -= dW_k
        self.W_v -= dW_v
        self.W_o -= dW_o

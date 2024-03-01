# Dylan Kenneth Eliot & GPT-4-Plugins ( Beta Edition )

"""

This basically lets one use the numpy & numba enough to optimize it without "optimizing it" for the root of all evil. This basically lets you run 3d calculations optimized to run entirely on the CPU and GPU where needed.

This basically lets you make use of the hardware minimally needed to use it. That's as it has been updated for
 as little as possible to prevent excessive optimization.

The way to make use of this copy is by updating your updated_main.py

Then run 

``docker run -it --rm -v $PWD/pybrain_networks/network\ pieces/:/app/ de3343/ai_mods_py3.10:numba_included sh -c "mv /string.py /usr/local/lib/python3.10/string.py; python3 
 /app/updated_main.py"``

to make use of it during hand modeling your version.


"""


from pybrain3.structure.modules.neuronlayer import NeuronLayer
from concurrent.futures import ThreadPoolExecutor
from numba import prange, jit
import numpy as np

@jit(parallel=True, nopython=True, nogil=True)
def manual_max(arr, axis=-1):
    if axis != -1:
        raise ValueError("This implementation only supports axis=-1.")
    
    if arr.ndim != 3:
        raise ValueError("This implementation only supports 3D arrays.")

    batch_size, seq_length, _ = arr.shape
    max_vals = np.empty((batch_size, seq_length), dtype=arr.dtype)

    for i in prange(batch_size):
        for j in prange(seq_length):
            max_val = arr[i, j, 0]
            for k in prange(arr.shape[2]):
                if arr[i, j, k] > max_val:
                    max_val = arr[i, j, k]
            max_vals[i, j] = max_val
    return max_vals.reshape(batch_size, seq_length, 1)


@jit(parallel=True, looplift=True)
def compute_attention_for_head(args):
    Q, K, V, dk = args
    K_transposed = np.transpose(K, axes=(0, 2, 1))
    batch_size, seq_length, _ = Q.shape
    matmul_qk = np.empty((batch_size, seq_length, seq_length), dtype=np.float64)
    for i in prange(batch_size):
        matmul_qk[i] = np.dot(Q[i], K_transposed[i])
    scaled_attention_logits = matmul_qk / np.sqrt(dk)
    # Use the manual_max function instead of np.max
    logits_max = manual_max(scaled_attention_logits, axis=-1)
    exp_logits = np.exp(scaled_attention_logits - logits_max)
    sum_exp_logits = np.sum(exp_logits, axis=-1).reshape(batch_size, seq_length, 1)
    attention_weights = exp_logits / sum_exp_logits
    output = np.empty_like(V)
    for i in prange(batch_size):
        output[i] = np.dot(attention_weights[i], V[i])
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
        args_list = [(Q[:, i, :, :], K[:, i, :, :], V[:, i, :, :], dk) for i in prange(Q.shape[1])]
        
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
    
    @jit(nopython=True, parallel=True, nogil=True)
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

        for i in prange(self.num_heads):
            dattention_output = dout_attention[:, :, i, :]
            dW_q_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))
            dW_k_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))
            dW_v_partial = np.dot(inbuf.reshape(-1, self.indim).T, dattention_output.reshape(-1, self.depth))

            dW_q += dW_q_partial
            dW_k += dW_k_partial
            dW_v += dW_v_partial

        self.W_q -= dW_q
        self.W_k -= dW_k
        self.W_v -= dW_v
        self.W_o -= dW_o
        dattention_output = np.concatenate(dattention_output, axis=-1)
        outbuf[:] = np.dot(dattention_output.reshape(batch_size, -1), self.W_o)

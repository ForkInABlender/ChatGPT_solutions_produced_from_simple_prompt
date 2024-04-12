# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""

Now there is a multi-dimensionally compatible encoding layer & 3d embeddign layer.

This is meant to run on the cpu only of their devices, meaning it is meant to be more or less cpu & memory based, instead of gpu based.
This should keep the scope to a managable size. Meaning that it should be usable in a model using 3-d attention head layers.




"""

from pybrain3.structure.modules import Module
import numpy as np

def softmax(x, axis=None):
    e_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
    return e_x / np.sum(e_x, axis=axis, keepdims=True)

class EmbeddingLayer3D(Module):
    def __init__(self, num_embeddings, embedding_dim):
        super(EmbeddingLayer3D, self).__init__()
        self.num_embeddings = num_embeddings
        self.embedding_dim = embedding_dim
        self.embeddings = np.random.randn(num_embeddings, embedding_dim) * 0.1
        # Additional array to store gradients for embeddings
        self.grad_embeddings = np.zeros_like(self.embeddings)

    def _forwardImplementation(self, inbuf, outbuf):
        assert np.all((inbuf >= 0) & (inbuf < self.num_embeddings)), "Indices out of bounds"
        embedded = self.embeddings[inbuf.astype(int)]
        outbuf[:] = embedded.flatten()

    def _backwardImplementation(self, outerr, inerr, inbuf):
        # Reshape outerr back to the shape of the output from forward pass
        reshaped_outerr = outerr.reshape(-1, self.embedding_dim)
        # Initialize gradient accumulation for embeddings
        self.grad_embeddings.fill(0)
        # Accumulate gradients for embeddings
        for i, idx in enumerate(inbuf.astype(int)):
            self.grad_embeddings[idx] += reshaped_outerr[i]
        # Pass gradient back to the input error (note: embeddings are normally end-layers, input error might not be used)
        inerr.fill(0)  # Typically, there is no backprop to indices, so this is just for compatibility

class Dim3EncodingLayer(Module):
    def __init__(self, features):
        super(GPTLikeEncodingLayer, self).__init__()
        self.features = features
        self.weights_q = np.random.randn(features, features) * 0.1
        self.weights_k = np.random.randn(features, features) * 0.1
        self.weights_v = np.random.randn(features, features) * 0.1

    def _reshape_input(self, inbuf):
        if inbuf.ndim == 1:
            return inbuf[None, None, :], (inbuf.shape[0],)
        elif inbuf.ndim == 2:
            return inbuf[None, :, :], (inbuf.shape[0], inbuf.shape[1])
        elif inbuf.ndim == 3:
            return inbuf, (inbuf.shape[0], inbuf.shape[1], inbuf.shape[2])
        else:
            raise ValueError("Input must be 1D, 2D, or 3D")

    def _reshape_output(self, outbuf, original_shape):
        if len(original_shape) == 1:
            return outbuf[0, 0, :]
        elif len(original_shape) == 2:
            return outbuf[0, :, :]
        elif len(original_shape) == 3:
            return outbuf
        else:
            raise ValueError("Output reshaping error")

    def _forwardImplementation(self, inbuf, outbuf):
        reshaped_input, original_shape = self._reshape_input(inbuf)
        query = np.einsum('ijk,kl->ijl', reshaped_input, self.weights_q)
        key = np.einsum('ijk,kl->ijl', reshaped_input, self.weights_k)
        value = np.einsum('ijk,kl->ijl', reshaped_input, self.weights_v)
        
        # Self-attention
        scale = self.features ** 0.5
        scores = np.matmul(query, key.transpose((0, 2, 1))) / scale
        weights = softmax(scores, axis=-1)
        attended = np.matmul(weights, value)
        outbuf[:] = self._reshape_output(attended, original_shape)

    def _backwardImplementation(self, outerr, inerr, inbuf):
        reshaped_input, original_shape = self._reshape_input(inbuf)
        query = np.einsum('ijk,kl->ijl', reshaped_input, self.weights_q)
        key = np.einsum('ijk,kl->ijl', reshaped_input, self.weights_k)
        value = np.einsum('ijk,kl->ijl', reshaped_input, self.weights_v)

        # Forward pass for reference
        scale = self.features ** 0.5
        scores = np.matmul(query, key.transpose((0, 2, 1))) / scale
        weights = softmax(scores, axis=-1)
        attended_output = np.matmul(weights, value)
        
        # Gradients of output w.r.t. weighted values
        d_out_d_attended = outerr.reshape(attended_output.shape)

        # Gradient w.r.t. weights used in weighted sum
        d_weights = np.matmul(d_out_d_attended, value.transpose((0, 2, 1)))

        # Gradient w.r.t. value
        d_value = np.matmul(weights.transpose((0, 2, 1)), d_out_d_attended)

        # Gradient w.r.t. scores (before softmax)
        d_scores = d_weights * (weights * (1 - weights))

        # Backprop through scale and transpose in scores computation
        d_query = np.matmul(d_scores / scale, key)
        d_key = np.matmul(d_scores.transpose((0, 2, 1)) / scale, query)

        # Backprop to input (reshaped_input) via query, key, value transformations
        d_reshaped_input_q = np.einsum('ijk,kl->ijl', d_query, self.weights_q.T)
        d_reshaped_input_k = np.einsum('ijk,kl->ijl', d_key, self.weights_k.T)
        d_reshaped_input_v = np.einsum('ijk,kl->ijl', d_value, self.weights_v.T)
        d_reshaped_input = d_reshaped_input_q + d_reshaped_input_k + d_reshaped_input_v

        # Gradient w.r.t. weights
        self.grad_weights_q = np.einsum('ijk,ijl->kl', reshaped_input, d_query)
        self.grad_weights_k = np.einsum('ijk,ijl->kl', reshaped_input, d_key)
        self.grad_weights_v = np.einsum('ijk,ijl->kl', reshaped_input, d_value)

        # Reshape gradient to match input shape for error propagation
        inerr[:] = self._reshape_output(d_reshaped_input, original_shape)

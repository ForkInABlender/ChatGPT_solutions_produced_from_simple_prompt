# Dylan Kenneth Eliot & GPT-4o ( alpha Edition )

"""

WARNING: TEST IN A AIR-GAPPED NETWORK OR OFFLINE IS RECOMMENDED UNLESS YOU KNOW WHAT YOU ARE DOING.

The given code could be refactored with { https://github.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/blob/2024_02/pybrain3_custom_layers/main__GPT_Model.py }
 for a more functional version.

This otherwise should be tested offline with unique datasets or otherwise real data, preferrably.

And practically, using tiktoken offline also would be useful. Now, what I am likely to do is use numpy, numba cuda simulator, docker, a google spreadsheet, and a flask server,
 Imitating GPT with pybrain3 should be a little bit of elbow grease.
"""

import numpy as np

def initialize_parameters(layer_dims):
    np.random.seed(3)
    parameters = {}
    L = len(layer_dims)

    for l in range(1, L):
        parameters['W' + str(l)] = np.random.randn(layer_dims[l], layer_dims[l-1]) * 0.01
        parameters['b' + str(l)] = np.zeros((layer_dims[l], 1))

    return parameters

def relu(Z):
    return np.maximum(0, Z)

def softmax(Z):
    expZ = np.exp(Z - np.max(Z))
    return expZ / expZ.sum(axis=0, keepdims=True)

def layer_normalization(X, epsilon=1e-6):
    mean = np.mean(X, axis=-1, keepdims=True)
    std = np.std(X, axis=-1, keepdims=True)
    return (X - mean) / (std + epsilon)

def linear_forward(A, W, b):
    return np.dot(W, A) + b

def self_attention(Q, K, V, mask=None):
    d_k = Q.shape[-1]
    scores = np.dot(Q, K.T) / np.sqrt(d_k)

    if mask is not None:
        scores = np.where(mask == 0, -1e9, scores)

    attention_weights = softmax(scores)
    return np.dot(attention_weights, V), attention_weights

def multi_head_attention(X, W_q, W_k, W_v, W_o, num_heads=8):
    d_model = X.shape[0]
    d_k = d_v = d_model // num_heads

    Q = np.dot(W_q, X).reshape((num_heads, d_k, -1))
    K = np.dot(W_k, X).reshape((num_heads, d_k, -1))
    V = np.dot(W_v, X).reshape((num_heads, d_v, -1))

    heads = []
    for i in range(num_heads):
        head, _ = self_attention(Q[i], K[i], V[i])
        heads.append(head)

    multi_head_output = np.concatenate(heads, axis=0)
    return np.dot(W_o, multi_head_output)

def feed_forward_network(X, W1, b1, W2, b2):
    return linear_forward(relu(linear_forward(X, W1, b1)), W2, b2)

def transformer_block(X, parameters, num_heads=8):
    W_q, W_k, W_v, W_o = parameters['W_q'], parameters['W_k'], parameters['W_v'], parameters['W_o']
    W1, b1, W2, b2 = parameters['W1'], parameters['b1'], parameters['W2'], parameters['b2']

    # Multi-Head Attention
    attn_output = multi_head_attention(X, W_q, W_k, W_v, W_o, num_heads)
    attn_output_norm = layer_normalization(attn_output + X)

    # Feed-Forward Network
    ff_output = feed_forward_network(attn_output_norm, W1, b1, W2, b2)
    ff_output_norm = layer_normalization(ff_output + attn_output_norm)

    return ff_output_norm

def positional_encoding(X):
    _, seq_len = X.shape
    pos = np.arange(seq_len)[:, np.newaxis]
    i = np.arange(X.shape[0])[np.newaxis, :]
    angle_rads = pos / np.power(10000, (2 * (i // 2)) / np.float32(X.shape[0]))

    angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])
    angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

    return X + angle_rads.T

def gpt_forward(X, parameters, num_layers=2, num_heads=8):
    X = positional_encoding(X)

    for _ in range(num_layers):
        X = transformer_block(X, parameters, num_heads)

    return softmax(X)

# Example Initialize parameters for the transformer block
layer_dims = [512, 2048]
parameters = {
    'W_q': np.random.randn(layer_dims[0], layer_dims[0]) * 0.01,
    'W_k': np.random.randn(layer_dims[0], layer_dims[0]) * 0.01,
    'W_v': np.random.randn(layer_dims[0], layer_dims[0]) * 0.01,
    'W_o': np.random.randn(layer_dims[0], layer_dims[0]) * 0.01,
    'W1': np.random.randn(layer_dims[1], layer_dims[0]) * 0.01,
    'b1': np.zeros((layer_dims[1], 1)),
    'W2': np.random.randn(layer_dims[0], layer_dims[1]) * 0.01,
    'b2': np.zeros((layer_dims[0], 1))
}

# Example usage
X = np.random.randn(512, 10)  # Input data with embedding dimension of 512 and sequence length of 10
output = gpt_forward(X, parameters, num_layers=2, num_heads=8)
print(output)

# Dylan Kenneth Eliot & GPT-4-Plugins ( Alpha Edition )

"""

This does the gpt-2 model definition without all the tensorflow boilerplate.

It should run just as dynamically as the one that openAI built. 


"""


import numpy as np

def default_hparams():
    return {
        'n_vocab': 0,
        'n_ctx': 1024,
        'n_embd': 768,
        'n_head': 12,
        'n_layer': 12,
    }

def shape_list(x):
    return list(x.shape)

def softmax(x, axis=-1):
    x = x - np.max(x, axis=axis, keepdims=True)
    ex = np.exp(x)
    return ex / np.sum(ex, axis=axis, keepdims=True)

def gelu(x):
    return 0.5 * x * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * np.power(x, 3))))

def norm(x, axis=-1, epsilon=1e-5):
    u = np.mean(x, axis=axis, keepdims=True)
    s = np.mean(np.square(x - u), axis=axis, keepdims=True)
    x = (x - u) * np.reciprocal(np.sqrt(s + epsilon))
    return x

def split_states(x, n):
    *start, m = x.shape
    return np.reshape(x, start + [n, m // n])

def merge_states(x):
    *start, a, b = x.shape
    return np.reshape(x, start + [a * b])

def conv1d(x, nf, w_init_stdev=0.02):
    *start, nx = x.shape
    w = np.random.normal(0, w_init_stdev, [1, nx, nf])
    b = np.zeros([nf])
    c = np.reshape(np.dot(np.reshape(x, [-1, nx]), np.reshape(w, [-1, nf])) + b, start + [nf])
    return c

def attention_mask(nd, ns, dtype=np.float32):
    i = np.arange(nd)[:, None]
    j = np.arange(ns)
    m = i >= j - ns + nd
    return m.astype(dtype)

def attn(x, n_state, hparams, past=None):
    assert x.ndim == 3
    assert n_state % hparams['n_head'] == 0

    def split_heads(x):
        return np.transpose(split_states(x, hparams['n_head']), [0, 2, 1, 3])

    def merge_heads(x):
        return merge_states(np.transpose(x, [0, 2, 1, 3]))

    def mask_attn_weights(w):
        _, _, nd, ns = w.shape
        b = attention_mask(nd, ns, dtype=w.dtype)
        b = np.reshape(b, [1, 1, nd, ns])
        w = w * b - np.cast[np.float32](1e10) * (1 - b)
        return w

    def multihead_attn(q, k, v):
        w = np.matmul(q, k.transpose([0, 1, 3, 2]))
        w = w * np.reciprocal(np.sqrt(v.shape[-1]))

        w = mask_attn_weights(w)
        w = softmax(w)
        a = np.matmul(w, v)
        return a

    c = conv1d(x, n_state * 3)
    q, k, v = np.split(c, 3, axis=2)
    q, k, v = map(split_heads, [q, k, v])
    if past is not None:
        pk, pv = np.split(past, 2, axis=1)
        k = np.concatenate([pk, k], axis=-2)
        v = np.concatenate([pv, v], axis=-2)
    a = multihead_attn(q, k, v)
    a = merge_heads(a)
    a = conv1d(a, n_state)
    return a

def mlp(x, n_state, hparams):
    nx = x.shape[-1]
    h = gelu(conv1d(x, n_state))
    h2 = conv1d(h, nx)
    return h2

def block(x, hparams, past=None):
    nx = x.shape[-1]
    a = attn(norm(x), nx, past=past, hparams=hparams)
    x = x + a
    m = mlp(norm(x), nx * 4, hparams=hparams)
    x = x + m
    return x

def model(hparams, X, past=None):
    # Assuming embedding matrices and positional encodings are initialized elsewhere
    wpe = np.random.normal(0, 0.01, [hparams['n_ctx'], hparams['n_embd']])  # Positional embeddings
    wte = np.random.normal(0, 0.02, [hparams['n_vocab'], hparams['n_embd']])  # Token embeddings
    
    # Embedding lookup for input tokens X and adding positional encodings
    past_length = 0 if past is None else past.shape[-2]
    position_encodings = np.arange(past_length, past_length + X.shape[1])
    h = wte[X, :] + wpe[position_encodings, :]
    
    # Transformer blocks
    presents = []
    for layer in range(hparams['n_layer']):
        h, present = block(h, hparams, past=past if past is not None else None)
        presents.append(present)
    
    # Final layer normalization (assuming norm is implemented correctly)
    h = norm(h)
    
    # Project back to vocabulary size for logits (assuming projection matrix is initialized elsewhere)
    projection_matrix = np.random.normal(0, 0.02, [hparams['n_embd'], hparams['n_vocab']])
    logits = np.dot(h, projection_matrix)
    
    # Compile results
    results = {
        'logits': logits,
        'present': np.stack(presents, axis=1)
    }
    
    return results

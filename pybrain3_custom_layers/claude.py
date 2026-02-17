#import sys
#sys.path.insert(0, '.')  # ensure string.py and pybrain3/ are on path

import scipy, numpy as np, types
"""_P = {
    'array':np.array,'asarray':np.asarray,'zeros':np.zeros,'ones':np.ones,
    'dot':np.dot,'outer':np.outer,'exp':np.exp,'tanh':np.tanh,'sqrt':np.sqrt,
    'log':np.log,'concatenate':np.concatenate,'reshape':np.reshape,
    'ravel':np.ravel,'size':np.size,'sum':np.sum,'mean':np.mean,'var':np.var,
    'where':np.where,'isinf':np.isinf,'isnan':np.isnan,'isscalar':np.isscalar,
    'ndarray':np.ndarray,'mat':np.matrix,'matrix':np.matrix,
    'randn':np.random.randn,'rand':np.random.rand,'random':np.random,
    'zeros_like':np.zeros_like,'empty':np.empty,'diag':np.diag,'eye':np.eye,
    'tile':np.tile,'arange':np.arange,'r_':np.r_,'c_':np.c_,'inf':np.inf,
    'sign':np.sign,'clip':np.clip,'multiply':np.multiply,'power':np.power,
    'pi':np.pi,'sin':np.sin,'cos':np.cos,'tan':np.tan,'floor':np.floor,
    'sort':np.sort,'argmax':np.argmax,'amin':np.amin,'amax':np.amax,
    'cov':np.cov,'asmatrix':np.asmatrix,'resize':np.resize,'arctanh':np.arctanh,
    'mgrid':np.mgrid,'log10':np.log10,'product':np.prod,
}
for k, v in _P.items():
    if not hasattr(scipy, k):
        setattr(scipy, k, v)
_w = types.ModuleType('scipy.weave')
sys.modules['scipy.weave'] = _w
setattr(scipy, 'weave', _w)
"""

from pybrain3.structure import SoftmaxLayer, LinearLayer, FeedForwardNetwork, FullConnection
from pybrain3.structure.modules.neuronlayer import NeuronLayer
from pybrain3.structure.modules.module import Module
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.datasets import SupervisedDataSet

# ─────────────────────────────────────────────────────────────
# GRADIENT CLIPPING
# ─────────────────────────────────────────────────────────────
CLIP = 1.0

def clip_grad(g):
    norm = np.linalg.norm(g)
    return g * (CLIP / norm) if norm > CLIP else g


# ─────────────────────────────────────────────────────────────
# CUSTOM LAYERS
# ─────────────────────────────────────────────────────────────

class FeedForwardLayer(Module):
    """Linear projection with bias. Uses Xavier init + gradient clipping."""
    def __init__(self, indim, outdim, name=None):
        super().__init__(indim, outdim, name=name)
        scale = np.sqrt(2.0 / (indim + outdim))
        self.weights = np.random.randn(indim, outdim) * scale
        self.bias    = np.zeros(outdim)
        self.lr      = 0.001

    def _forwardImplementation(self, inbuf, outbuf):
        outbuf[:] = np.clip(inbuf @ self.weights + self.bias, -10, 10)

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        outerr = clip_grad(outerr)
        dW = inbuf[:, np.newaxis] @ outerr[np.newaxis, :]
        self.weights -= self.lr * clip_grad(dW)
        self.bias    -= self.lr * outerr
        inerr[:]      = clip_grad(outerr @ self.weights.T)


class EmbeddingLayer(Module):
    """Token index (one-hot input) → dense embedding vector."""
    def __init__(self, vocab_size, embedding_dim, name=None):
        super().__init__(vocab_size, embedding_dim, name=name)
        self.embeddings = np.random.randn(vocab_size, embedding_dim) * 0.01
        self.token_idx  = 0
        self.lr         = 0.001

    def _forwardImplementation(self, inbuf, outbuf):
        self.token_idx = int(np.argmax(inbuf))
        outbuf[:] = self.embeddings[min(self.token_idx, len(self.embeddings) - 1)]

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        grad = np.zeros_like(self.embeddings)
        grad[self.token_idx] = clip_grad(outerr)
        self.embeddings -= self.lr * grad
        inerr[:] = 0


class GeLULayer(NeuronLayer):
    """GELU activation — exact GPT formula."""
    def __init__(self, dim, name=None):
        super().__init__(dim, name=name)

    def _forwardImplementation(self, inbuf, outbuf):
        x = np.clip(inbuf, -10, 10)
        outbuf[:] = 0.5 * x * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x ** 3)))

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        x   = np.clip(inbuf, -10, 10)
        cdf = 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x ** 3)))
        arg = np.clip(np.sqrt(2 / np.pi) * (x + 0.044715 * x ** 3), -30, 30)
        pd  = np.sqrt(2 / np.pi) * (1 + 0.134145 * x ** 2) * (1 / np.cosh(arg)) ** 2
        inerr[:] = clip_grad(outerr * (cdf + x * pd))


def softmax(x, axis=-1):
    e = np.exp(np.clip(x - np.max(x, axis=axis, keepdims=True), -30, 0))
    return e / (e.sum(axis=axis, keepdims=True) + 1e-9)


class MultiHeadSelfAttention(NeuronLayer):
    """Scaled dot-product multi-head self-attention."""
    def __init__(self, dim, num_heads, name=None):
        super().__init__(dim, name=name)
        self.num_heads = num_heads
        scale = np.sqrt(2.0 / (dim * 2))
        self.W_q = np.random.randn(dim, dim) * scale
        self.W_k = np.random.randn(dim, dim) * scale
        self.W_v = np.random.randn(dim, dim) * scale
        self.W_o = np.random.randn(dim, dim) * scale
        self.lr  = 0.001

    def scaled_dot_product_attention(self, Q, K, V):
        scores = np.dot(Q, K.T) / np.sqrt(max(Q.shape[-1], 1))
        return np.dot(softmax(scores, axis=-1), V)

    def _forwardImplementation(self, inbuf, outbuf):
        x  = inbuf[np.newaxis, :] if inbuf.ndim == 1 else inbuf
        Qs = np.split(np.dot(x, self.W_q), self.num_heads, axis=1)
        Ks = np.split(np.dot(x, self.W_k), self.num_heads, axis=1)
        Vs = np.split(np.dot(x, self.W_v), self.num_heads, axis=1)
        heads = [self.scaled_dot_product_attention(Qs[i], Ks[i], Vs[i])
                 for i in range(self.num_heads)]
        out = np.dot(np.concatenate(heads, axis=1), self.W_o)
        outbuf[:] = clip_grad(out.ravel()[:outbuf.shape[0]])

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        inerr[:] = clip_grad(outerr @ self.W_o.T[:outerr.shape[0], :inerr.shape[0]])


class LayerNorm(NeuronLayer):
    """Layer normalisation with learnable gamma/beta."""
    def __init__(self, dim, eps=1e-6, name=None):
        super().__init__(dim, name=name)
        self.gamma = np.ones(dim)
        self.beta  = np.zeros(dim)
        self.eps   = eps

    def _forwardImplementation(self, inbuf, outbuf):
        mu  = np.mean(inbuf)
        sd  = np.std(inbuf) + self.eps
        outbuf[:] = self.gamma * (inbuf - mu) / sd + self.beta

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        N  = inbuf.size
        sd = np.std(inbuf) + self.eps
        inerr[:] = clip_grad(
            self.gamma / (N * sd) * (
                N * outerr - np.sum(outerr) -
                (inbuf - np.mean(inbuf)) / sd * np.sum(outerr * (inbuf - np.mean(inbuf)))
            )
        )


# ─────────────────────────────────────────────────────────────
# TOKENISER  (character-level, printable ASCII)
# ─────────────────────────────────────────────────────────────
import string as _str
VOCAB_CHARS = list(_str.printable)
VOCAB_SIZE  = len(VOCAB_CHARS)          # 100
char2idx    = {c: i for i, c in enumerate(VOCAB_CHARS)}
idx2char    = {i: c for c, i in char2idx.items()}

def tokenise(text):
    return [char2idx.get(c, 0) for c in text]

def detokenise(ids):
    return ''.join(idx2char.get(i, '?') for i in ids)

def one_hot(i, n=VOCAB_SIZE):
    v = np.zeros(n)
    v[i] = 1.0
    return v


# ─────────────────────────────────────────────────────────────
# ARCHITECTURE
#   Real Claude:  d_model~8192, ~80 layers, ~32 heads, 100K vocab
#   This model:   d_model=64,   4 layers,  4 heads,   100 vocab
#   Structure is identical — only scale differs.
# ─────────────────────────────────────────────────────────────
D      = 64
HEADS  = 4    # D / HEADS = 16 per head
BLOCKS = 4
FFN_D  = 64


def build_net(d=D, heads=HEADS, blocks=BLOCKS, ffn_d=FFN_D, vocab=VOCAB_SIZE):
    net  = FeedForwardNetwork()
    inp  = LinearLayer(vocab, name='token_input')
    emb  = EmbeddingLayer(vocab, d, name='embedding')
    mhsa = MultiHeadSelfAttention(d, heads, name='mhsa')

    net.addInputModule(inp)
    net.addModule(emb);  net.addModule(mhsa)
    net.addConnection(FullConnection(inp, emb))
    net.addConnection(FullConnection(emb, mhsa))

    prev = mhsa
    for i in range(blocks):
        ln1  = LayerNorm(d,        name=f'ln1_{i}')
        ff1  = FeedForwardLayer(d, ffn_d, name=f'ff1_{i}')
        gelu = GeLULayer(ffn_d,    name=f'gelu_{i}')
        ff2  = FeedForwardLayer(ffn_d, d, name=f'ff2_{i}')
        ln2  = LayerNorm(d,        name=f'ln2_{i}')
        for m in [ln1, ff1, gelu, ff2, ln2]:
            net.addModule(m)
        net.addConnection(FullConnection(prev, ln1))
        net.addConnection(FullConnection(ln1,  ff1))
        net.addConnection(FullConnection(ff1,  gelu))
        net.addConnection(FullConnection(gelu, ff2))
        net.addConnection(FullConnection(ff2,  ln2))
        prev = ln2

    out = SoftmaxLayer(vocab, name='output')
    net.addOutputModule(out)
    net.addConnection(FullConnection(prev, out))
    net.sortModules()
    return net


# ─────────────────────────────────────────────────────────────
# TRAINING DATA
# Replace CORPUS with your own text to train on anything.
# ─────────────────────────────────────────────────────────────
CORPUS = """What is 2+2? 2+2 equals 4. What is the capital of France? The capital of France is Paris. What is Python? Python is a high-level programming language. What is gravity? Gravity is a force attracting objects with mass.""".strip()


def make_dataset(text):
    tokens = tokenise(text)
    ds = SupervisedDataSet(VOCAB_SIZE, VOCAB_SIZE)
    for i in range(len(tokens) - 1):
        ds.addSample(one_hot(tokens[i]), one_hot(tokens[i + 1]))
    return ds, tokens


# ─────────────────────────────────────────────────────────────
# BUILD + TRAIN
# ─────────────────────────────────────────────────────────────
net = build_net()
ds, tokens = make_dataset(CORPUS)

print("╔══════════════════════════════════════════════════════════╗")
print("║       Claude via pybrain3 — full pipeline                ║")
print("╚══════════════════════════════════════════════════════════╝")
print(f"  vocab={VOCAB_SIZE}  d={D}  heads={HEADS}  blocks={BLOCKS}")
print(f"  params={len(net.params):,}  corpus={len(CORPUS)} chars  samples={len(ds)}")

trainer = BackpropTrainer(net, ds, learningrate=0.002, momentum=0.9, verbose=False)

import warnings
warnings.filterwarnings('ignore')

EPOCHS = 100
print(f"\n  Training {EPOCHS} epochs (next-char prediction)...")
errors = []
for epoch in range(EPOCHS):
    err = trainer.train()
    errors.append(err)
    if (epoch + 1) % 20 == 0:
        finite = [e for e in errors[-20:] if np.isfinite(e)]
        avg = np.mean(finite) if finite else float('nan')
        print(f"    epoch {epoch + 1:>3}  avg_error={avg:.6f}")

finite_errors = [e for e in errors if np.isfinite(e)]
print(f"\n  Finite error samples : {len(finite_errors)}/{len(errors)}")
if finite_errors:
    print(f"  Error start  : {finite_errors[0]:.6f}")
    print(f"  Error end    : {finite_errors[-1]:.6f}")
    print(f"  Improved     : {finite_errors[-1] < finite_errors[0]}")


# ─────────────────────────────────────────────────────────────
# INFERENCE — temperature-sampled next-char generation
# ─────────────────────────────────────────────────────────────
def generate(net, seed, n=50, temperature=0.8):
    tokens = tokenise(seed)
    out_chars = list(seed)
    cur = tokens[-1]
    for _ in range(n):
        probs = net.activate(one_hot(cur))
        probs = np.clip(probs, 1e-9, None)
        probs = probs ** (1.0 / temperature)
        probs /= probs.sum()
        if np.any(np.isnan(probs)):
            probs = np.ones(VOCAB_SIZE) / VOCAB_SIZE
        nxt = int(np.random.choice(len(probs), p=probs))
        out_chars.append(idx2char.get(nxt, '?'))
        cur = nxt
    return ''.join(out_chars)


print(f"\n  Generation (temperature=0.8):")
for seed in ["What is", "Python is", "The capital"]:
    result = generate(net, seed, n=40, temperature=0.8)
    print(f"    [{seed!r}] → '{result}'")


# ─────────────────────────────────────────────────────────────
# SAVE SNAPSHOT
# ─────────────────────────────────────────────────────────────
from pybrain3.tools.xml.networkwriter import NetworkWriter
from pybrain3.tools.xml.networkreader import NetworkReader

save_path = 'claude_pybrain3_snapshot.xml'
NetworkWriter.writeToFile(net, save_path)
print(f"\n  Saved → {save_path}")
print(f"  Reload with: net = NetworkReader.readFrom('{save_path}')")

# Dylan Kenneth Eliot


"""
Knowledge distillation: n-gram teacher -> LSTM student.

Phase 1: Train LSTM on soft targets from n-gram (distillation)
Phase 2: Fine-tune LSTM on hard corpus targets (ground truth)

Imports shared primitives from vortex_lm.py.
"""

import numpy as np
from collections import defaultdict, Counter
from sklearn.cluster import SpectralClustering
from pybrain3.structure import RecurrentNetwork, LSTMLayer, LinearLayer, SoftmaxLayer, FullConnection
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.datasets import SequentialDataSet

# ── Copy shared primitives (avoid import side-effects from vortex_lm) ──────────

DOUBLING_SEQ = [1, 2, 4, 8, 7, 5]
DOUBLING_POS = {v: i for i, v in enumerate(DOUBLING_SEQ)}

def digital_root(x):
    return (x % 9) or 9

def vortex_pos(x):
    return DOUBLING_POS.get(digital_root(x), -1)

def family_group(x):
    dr = digital_root(x)
    if dr in {1, 4, 7}: return 0
    if dr in {2, 5, 8}: return 1
    return 2

def build_vocab(corpus):
    """
    Contiguous char indices for net I/O.
    Vortex cluster membership stored separately as features only —
    not baked into token IDs, so no single output neuron is privileged.
    """
    chars = sorted(set(corpus))
    n = len(chars)
    c2i_raw = {c: i for i, c in enumerate(chars)}

    W = np.zeros((n, n))
    for k in range(len(corpus)):
        for j in range(k + 1, min(k + 3, len(corpus))):
            a, b = c2i_raw[corpus[k]], c2i_raw[corpus[j]]
            W[a, b] += 1; W[b, a] += 1

    n_clusters = min(9, n)
    labels = SpectralClustering(
        n_clusters=n_clusters, affinity='precomputed', random_state=42
    ).fit_predict(W + 1e-6 * np.eye(n))

    char2id = {c: i for i, c in enumerate(chars)}   # char -> 0..n-1
    id2char = {i: c for c, i in char2id.items()}
    char2cluster = {chars[i]: int(labels[i]) for i in range(n)}
    id2idx = {i: i for i in range(n)}   # identity: index == id
    idx2id = id2idx
    return char2id, id2char, id2idx, idx2id, char2cluster

def encode_char(vid, id2idx, vocab_size, pos, char2cluster=None, id2char=None):
    vec = np.zeros(vocab_size + 3)
    vec[id2idx[vid]] = 1.0
    # vortex features from cluster label
    ch  = id2char[vid] if id2char else None
    lbl = char2cluster[ch] if (char2cluster and ch) else 0
    vec[vocab_size]     = (lbl % 6) / 5.0          # doubling position proxy
    vec[vocab_size + 1] = (lbl // 3) / 2.0         # family group proxy
    vec[vocab_size + 2] = digital_root(pos) / 9.0  # positional dr
    return vec

# ── N-gram teacher ─────────────────────────────────────────────────────────────

class NGramTeacher:
    def __init__(self, n=4):
        self.n = n
        self.counts = defaultdict(Counter)

    def train(self, text):
        for i in range(len(text) - self.n):
            self.counts[text[i:i + self.n]][text[i + self.n]] += 1

    def soft_target(self, ctx, vocab, c2i, vocab_size):
        """Return soft probability vector over full vocab."""
        key = ctx[-(self.n):]
        while key and key not in self.counts:
            key = key[1:]
        y = np.zeros(vocab_size)
        if not key:
            return None
        c = self.counts[key]
        total = sum(c.values())
        for ch, cnt in c.items():
            if ch in c2i:
                y[c2i[ch]] = cnt / total
        return y

# ── Datasets ───────────────────────────────────────────────────────────────────

def distill_dataset(corpus, ngram, char2id, id2idx, vocab_size, char2cluster, id2char, chunk=15):
    ids  = [char2id[c] for c in corpus if c in char2id]
    text = [c for c in corpus if c in char2id]
    feat = vocab_size + 3
    ds   = SequentialDataSet(feat, vocab_size)
    c2i  = {ch: i for i, ch in enumerate(sorted(char2id.keys()))}
    for start in range(0, len(ids) - 1, chunk):
        seg_ids  = ids[start:start + chunk + 1]
        seg_text = text[start:start + chunk + 1]
        ds.newSequence()
        for i in range(len(seg_ids) - 1):
            x   = encode_char(seg_ids[i], id2idx, vocab_size, start + i + 1, char2cluster, id2char)
            ctx = ''.join(seg_text[max(0, i - ngram.n + 1):i + 1])
            y   = ngram.soft_target(ctx, sorted(char2id.keys()), c2i, vocab_size)
            if y is None:
                y = np.zeros(vocab_size); y[id2idx[seg_ids[i + 1]]] = 1.0
            ds.addSample(x, y)
    return ds

def hard_dataset(corpus, char2id, id2idx, vocab_size, char2cluster, id2char, chunk=15):
    ids  = [char2id[c] for c in corpus if c in char2id]
    feat = vocab_size + 3
    ds   = SequentialDataSet(feat, vocab_size)
    for start in range(0, len(ids) - 1, chunk):
        seg = ids[start:start + chunk + 1]
        ds.newSequence()
        for i in range(len(seg) - 1):
            x = encode_char(seg[i], id2idx, vocab_size, start + i + 1, char2cluster, id2char)
            y = np.zeros(vocab_size); y[id2idx[seg[i + 1]]] = 1.0
            ds.addSample(x, y)
    return ds

# ── RNN ────────────────────────────────────────────────────────────────────────

def build_rnn(input_size, hidden_size, output_size):
    net = RecurrentNetwork()
    net.addInputModule(LinearLayer(input_size,  name='in'))
    net.addModule(LSTMLayer(hidden_size,         name='lstm'))
    net.addOutputModule(SoftmaxLayer(output_size, name='out'))
    net.addConnection(FullConnection(net['in'],   net['lstm']))
    net.addConnection(FullConnection(net['lstm'], net['out']))
    net.addRecurrentConnection(FullConnection(net['lstm'], net['lstm'], name='rec'))
    net.sortModules()
    return net

def neural_generate(net, seed, char2id, id2char, id2idx, idx2id, vocab_size, char2cluster, length=200, temp=0.8):
    ids = [char2id[c] for c in seed if c in char2id]
    result = list(seed)
    net.reset()
    for pos, vid in enumerate(ids):
        net.activate(encode_char(vid, id2idx, vocab_size, pos + 1, char2cluster, id2char))
    pos, last_vid = len(ids), ids[-1]
    for _ in range(length):
        probs = np.array(net.activate(encode_char(last_vid, id2idx, vocab_size, pos, char2cluster, id2char)), dtype=np.float64)
        log_p = np.log(np.clip(probs, 1e-9, 1.0)) / temp
        log_p -= log_p.max()
        probs  = np.exp(log_p); probs /= probs.sum()
        last_vid = idx2id[np.random.choice(len(probs), p=probs)]
        result.append(id2char[last_vid])
        pos += 1
    return ''.join(result)

# ── Corpus ─────────────────────────────────────────────────────────────────────

CORPUS = """
the cat sat on the mat the cat ate the rat
a man a plan a canal panama
to be or not to be that is the question
all that glitters is not gold
the quick brown fox jumps over the lazy dog
she sells seashells by the seashore
how much wood would a woodchuck chuck
peter piper picked a peck of pickled peppers
i think therefore i am
the only way to do great work is to love what you do
in the beginning was the word and the word was with god
to infinity and beyond
elementary my dear watson
may the force be with you
the truth is out there
with great power comes great responsibility
it does not matter how slowly you go as long as you do not stop
the journey of a thousand miles begins with one step
you miss one hundred percent of the shots you do not take
whether you think you can or you think you cannot you are right
language is the road map of a culture
a language that does not affect the way you think is not worth knowing
the limits of my language mean the limits of my world
one language sets you in a corridor for life two languages open every door
to have another language is to possess a second soul
knowledge is power and power is knowledge
the mind is not a vessel to be filled but a fire to be kindled
we are what we repeatedly do excellence then is not an act but a habit
the unexamined life is not worth living
imagination is more important than knowledge
in the middle of difficulty lies opportunity
life is what happens when you are busy making other plans
the best way to predict the future is to create it
do not go where the path may lead go instead where there is no path
you cannot step into the same river twice
the only true wisdom is in knowing you know nothing
an unexamined life is not worth living
the measure of intelligence is the ability to change
logic will get you from a to b imagination will take you everywhere
""".strip()

# ── Run ────────────────────────────────────────────────────────────────────────

print("Building vocab...")
char2id, id2char, id2idx, idx2id, char2cluster = build_vocab(CORPUS)
vocab_size = len(char2id)
feat_size  = vocab_size + 3
print(f"Vocab size: {vocab_size}  Feat size: {feat_size}")

print("Training n-gram teacher...")
ngram = NGramTeacher(n=4)
ngram.train(CORPUS)

print("Building LSTM...")
net = build_rnn(feat_size, hidden_size=64, output_size=vocab_size)

epochs_=200

# ── Phase 1: distillation ──────────────────────────────────────────────────────
print(f"\nPhase 1: distill from n-gram teacher ({epochs_} epochs)...")
ds_soft = distill_dataset(CORPUS, ngram, char2id, id2idx, vocab_size, char2cluster, id2char)
trainer  = BackpropTrainer(net, ds_soft, learningrate=0.01, momentum=0.9)
for epoch in range(1, epochs_):
    err = trainer.train()
    if 0 == 0:
        print(f"  epoch {epoch:3d}  error={err:.6f}")

# ── Phase 2: blended fine-tune ─────────────────────────────────────────────────
print(f"\nPhase 2: fine-tune on blended targets ({epochs_} epochs)...")

def blended_dataset(corpus, ngram, char2id, id2idx, vocab_size, char2cluster, id2char, chunk=15, alpha=0.3):
    ids  = [char2id[c] for c in corpus if c in char2id]
    text = [c for c in corpus if c in char2id]
    feat = vocab_size + 3
    ds   = SequentialDataSet(feat, vocab_size)
    c2i  = {ch: i for i, ch in enumerate(sorted(char2id.keys()))}
    for start in range(0, len(ids) - 1, chunk):
        seg_ids  = ids[start:start + chunk + 1]
        seg_text = text[start:start + chunk + 1]
        ds.newSequence()
        for i in range(len(seg_ids) - 1):
            x    = encode_char(seg_ids[i], id2idx, vocab_size, start + i + 1, char2cluster, id2char)
            ctx  = ''.join(seg_text[max(0, i - ngram.n + 1):i + 1])
            soft = ngram.soft_target(ctx, sorted(char2id.keys()), c2i, vocab_size)
            hard = np.zeros(vocab_size); hard[id2idx[seg_ids[i + 1]]] = 1.0
            y    = alpha * hard + (1 - alpha) * (soft if soft is not None else hard)
            ds.addSample(x, y)
    return ds

ds_blend = blended_dataset(CORPUS, ngram, char2id, id2idx, vocab_size, char2cluster, id2char)
trainer   = BackpropTrainer(net, ds_blend, learningrate=0.0005, momentum=0.9)
for epoch in range(1, epochs_):
    err = trainer.train()
    if 0 == 0:
        print(f"  epoch {epoch:3d}  error={err:.6f}")

# ── Generate ───────────────────────────────────────────────────────────────────
seeds = ["the ", "language ", "to be "]
print("\n── Distilled LSTM+Vortex responses ──────────────────────────────")
for seed in seeds:
    print(f"\n[seed='{seed}']")
    print(neural_generate(net, seed, char2id, id2char, id2idx, idx2id,
                          vocab_size, char2cluster, length=150, temp=0.6))

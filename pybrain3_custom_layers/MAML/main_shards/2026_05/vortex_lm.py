"""
Vortex-augmented character-level LM.
- RecurrentNetwork with LSTMLayer + SoftmaxLayer (pybrain3)
- Vortex co-occurrence vocab alignment
- N-gram fallback generator for coherent output when neural output is weak
"""

import numpy as np
from collections import defaultdict, Counter
from sklearn.cluster import SpectralClustering
from pybrain3.structure import RecurrentNetwork, LSTMLayer, LinearLayer, SoftmaxLayer, FullConnection
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.datasets import SequentialDataSet

# ── Vortex primitives ──────────────────────────────────────────────────────────

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

# ── Vocab ──────────────────────────────────────────────────────────────────────

def build_vocab(corpus):
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

    cluster_weight = defaultdict(float)
    for i, lbl in enumerate(labels):
        cluster_weight[lbl] += W[i].sum()
    sorted_clusters = sorted(cluster_weight, key=cluster_weight.get, reverse=True)
    target_drs = DOUBLING_SEQ + [3, 6, 9]
    cluster_to_dr = {c: target_drs[i] for i, c in enumerate(sorted_clusters)}

    def ids_with_dr(dr, count):
        result, n = [], 1
        while len(result) < count:
            if digital_root(n) == dr: result.append(n)
            n += 1
        return result

    cluster_chars = defaultdict(list)
    for i, lbl in enumerate(labels):
        cluster_chars[lbl].append(chars[i])

    char2id, id2char = {}, {}
    for lbl, clist in cluster_chars.items():
        dr = cluster_to_dr[lbl]
        ids = ids_with_dr(dr, len(clist))
        for ch, vid in zip(clist, ids):
            char2id[ch] = vid
            id2char[vid] = ch

    all_ids = sorted(id2char.keys())
    id2idx  = {vid: i for i, vid in enumerate(all_ids)}
    idx2id  = {i: vid for vid, i in id2idx.items()}
    return char2id, id2char, id2idx, idx2id

# ── Encoding ───────────────────────────────────────────────────────────────────

def encode_char(vid, id2idx, vocab_size, pos):
    vec = np.zeros(vocab_size + 3)
    vec[id2idx[vid]] = 1.0
    vec[vocab_size]     = (vortex_pos(vid) + 1) / 6.0
    vec[vocab_size + 1] = family_group(vid) / 2.0
    vec[vocab_size + 2] = digital_root(pos) / 9.0
    return vec

# ── N-gram fallback ────────────────────────────────────────────────────────────

class NGramModel:
    """
    N-gram model with a pybrain3 SoftmaxLayer re-scorer.
    - counts table handles backoff lookup (unchanged)
    - a FeedForwardNetwork with SoftmaxLayer re-scores the full vocab
      distribution at each step using the n-gram context as input
    - generate samples from the re-scored distribution
    """
    def __init__(self, n=4):
        self.n = n
        self.counts = defaultdict(Counter)
        self.net = None
        self._vocab = None

    def train(self, text):
        for i in range(len(text) - self.n):
            ctx = text[i:i + self.n]
            nxt = text[i + self.n]
            self.counts[ctx][nxt] += 1

    def build_scorer(self, vocab):
        """Build and train a pybrain3 net to re-score n-gram distributions."""
        from pybrain3.structure import FeedForwardNetwork, LinearLayer, SoftmaxLayer, FullConnection
        from pybrain3.supervised.trainers import BackpropTrainer
        from pybrain3.datasets import SupervisedDataSet

        self._vocab = vocab
        v = len(vocab)
        c2i = {c: i for i, c in enumerate(vocab)}

        # input: one-hot of each char in context (n * v)
        # output: one-hot of next char (v) — trained with softmax
        inp_size = self.n * v
        net = FeedForwardNetwork()
        net.addInputModule(LinearLayer(inp_size, name='in'))
        net.addOutputModule(SoftmaxLayer(v, name='out'))
        net.addConnection(FullConnection(net['in'], net['out']))
        net.sortModules()

        ds = SupervisedDataSet(inp_size, v)
        for ctx, nexts in self.counts.items():
            if len(ctx) != self.n:
                continue
            x = np.zeros(inp_size)
            for j, c in enumerate(ctx):
                if c in c2i:
                    x[j * v + c2i[c]] = 1.0
            total = sum(nexts.values())
            y = np.zeros(v)
            for c, cnt in nexts.items():
                if c in c2i:
                    y[c2i[c]] = cnt / total
            ds.addSample(x, y)

        trainer = BackpropTrainer(net, ds, learningrate=0.05, momentum=0.9)
        for _ in range(30):
            trainer.train()

        self.net = net
        self._c2i = c2i
        self._i2c = {i: c for c, i in c2i.items()}

    def _rescore(self, ctx, raw_chars, raw_probs):
        """Re-score via pybrain3 softmax net if available, else return raw."""
        if self.net is None:
            return raw_chars, raw_probs
        v = len(self._vocab)
        x = np.zeros(self.n * v)
        padded = ctx[-(self.n):].rjust(self.n)
        for j, c in enumerate(padded):
            if c in self._c2i:
                x[j * v + self._c2i[c]] = 1.0
        scores = np.array(self.net.activate(x))
        # blend: 70% neural re-score, 30% raw n-gram
        raw_full = np.zeros(v)
        for c, p in zip(raw_chars, raw_probs):
            if c in self._c2i:
                raw_full[self._c2i[c]] = p
        blended = 0.3 * scores + 0.7 * raw_full
        blended /= blended.sum()
        return list(self._i2c.values()), list(blended)

    def predict(self, ctx):
        key = ctx[-(self.n):]
        while key and key not in self.counts:
            key = key[1:]
        if not key:
            return None
        c = self.counts[key]
        total = sum(c.values())
        chars = list(c.keys())
        probs = [v / total for v in c.values()]
        # only rescore if net is trained; keep chars restricted to n-gram candidates
        if self.net is None:
            return chars, probs
        v = len(self._vocab)
        x = np.zeros(self.n * v)
        padded = ctx[-(self.n):].rjust(self.n)
        for j, ch in enumerate(padded):
            if ch in self._c2i:
                x[j * v + self._c2i[ch]] = 1.0
        scores = np.array(self.net.activate(x))
        raw_full = np.zeros(v)
        for ch, p in zip(chars, probs):
            if ch in self._c2i:
                raw_full[self._c2i[ch]] = p
        blended = 0.3 * scores + 0.7 * raw_full
        # restrict back to n-gram candidates only to preserve coherence
        candidate_idx = [self._c2i[ch] for ch in chars if ch in self._c2i]
        restricted = np.zeros(v)
        for i in candidate_idx:
            restricted[i] = blended[i]
        total_r = restricted.sum()
        if total_r == 0:
            return chars, probs
        return [self._i2c[i] for i in candidate_idx], [restricted[i] / total_r for i in candidate_idx]

    def generate(self, seed, length=200, temp=0.8):
        result = list(seed)
        for _ in range(length):
            out = self.predict(''.join(result))
            if out is None:
                break
            chars, probs = out
            probs = np.array(probs) ** (1.0 / temp)
            probs /= probs.sum()
            result.append(np.random.choice(chars, p=probs))
        return ''.join(result)

# ── Neural net ─────────────────────────────────────────────────────────────────

def build_rnn(input_size, hidden_size, output_size):
    net = RecurrentNetwork()
    net.addInputModule(LinearLayer(input_size, name='in'))
    net.addModule(LSTMLayer(hidden_size, name='lstm'))
    net.addOutputModule(SoftmaxLayer(output_size, name='out'))
    net.addConnection(FullConnection(net['in'],   net['lstm']))
    net.addConnection(FullConnection(net['lstm'], net['out']))
    net.addRecurrentConnection(FullConnection(net['lstm'], net['lstm'], name='rec'))
    net.sortModules()
    return net

def build_dataset(text, char2id, id2idx, vocab_size, chunk=30):
    """Split into short chunks so BPTT doesn't explode."""
    ids = [char2id[c] for c in text if c in char2id]
    feat = vocab_size + 3
    ds = SequentialDataSet(feat, vocab_size)
    for start in range(0, len(ids) - 1, chunk):
        seg = ids[start:start + chunk + 1]
        ds.newSequence()
        for i in range(len(seg) - 1):
            x = encode_char(seg[i], id2idx, vocab_size, start + i + 1)
            y = np.zeros(vocab_size)
            y[id2idx[seg[i + 1]]] = 1.0
            ds.addSample(x, y)
    return ds

def neural_generate(net, seed, char2id, id2char, id2idx, idx2id, vocab_size, length=200, temp=0.8):
    ids = [char2id[c] for c in seed if c in char2id]
    result = list(seed)
    net.reset()

    # prime with seed
    for pos, vid in enumerate(ids):
        x = encode_char(vid, id2idx, vocab_size, pos + 1)
        raw = net.activate(x)
    
    pos = len(ids)
    last_vid = ids[-1]
    for _ in range(length):
        x = encode_char(last_vid, id2idx, vocab_size, pos)
        probs = np.array(net.activate(x), dtype=np.float64)
        # temperature scaling on log probs then renorm (claude.py style)
        log_p = np.log(np.clip(probs, 1e-9, 1.0)) / temp
        log_p -= log_p.max()
        probs = np.exp(log_p)
        probs /= probs.sum()
        idx = np.random.choice(len(probs), p=probs)
        last_vid = idx2id[idx]
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
char2id, id2char, id2idx, idx2id = build_vocab(CORPUS)
vocab_size = len(char2id)
feat_size  = vocab_size + 3
print(f"Vocab size: {vocab_size}  Feat size: {feat_size}")

print("Training n-gram model...")
ngram = NGramModel(n=4)
ngram.train(CORPUS)
print("Building n-gram pybrain3 re-scorer...")
ngram.build_scorer(sorted(set(CORPUS)))

print("Building RNN dataset...")
ds = build_dataset(CORPUS, char2id, id2idx, vocab_size)

print("Building recurrent network...")
net = build_rnn(feat_size, hidden_size=64, output_size=vocab_size)

print("Training RNN (30 epochs)...")
trainer = BackpropTrainer(net, ds, learningrate=0.001, momentum=0.9)
for epoch in range(1, 31):
    err = trainer.train()
    if epoch % 5 == 0:
        print(f"  epoch {epoch:3d}  error={err:.6f}")

seeds = ["the ", "language ", "to be "]

print("\n── N-gram responses ─────────────────────────────────────────────")
for seed in seeds:
    print(f"\n[seed='{seed}']")
    print(ngram.generate(seed, length=150, temp=0.8))

print("\n── Neural (LSTM+Vortex) responses ───────────────────────────────")
for seed in seeds:
    print(f"\n[seed='{seed}']")
    print(neural_generate(net, seed, char2id, id2char, id2idx, idx2id,
                          vocab_size, length=150, temp=0.8))

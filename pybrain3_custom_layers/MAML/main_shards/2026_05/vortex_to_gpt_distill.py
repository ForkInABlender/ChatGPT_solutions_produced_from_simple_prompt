# Dylan Kenneth Eliot
"""
Distill from vortex_distill.py LSTM teacher -> main_gpt3_3d.py GPT student.

Uses n-gram soft targets (matching vortex_distill.py's teacher pattern) over
GPT's 50257 token space. Vortex math encodes positional weight on the input
one-hot — no extra dimensions, 50257 I/O fully preserved.
"""

import numpy as np
from collections import defaultdict, Counter
from pybrain3.tools.xml import NetworkReader, NetworkWriter
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer

from vortex_distill import (
    neural_generate, CORPUS,
    char2id, id2char, id2idx, idx2id, char2cluster, vocab_size,
    digital_root, family_group, vortex_pos,
    NGramTeacher
)

VOCAB_SIZE = 50257

# ── load vortex LSTM teacher ───────────────────────────────────────────────────
lstm_teacher = NetworkReader.readFrom("vortex_model.xml")
print("Loaded vortex LSTM teacher.")

# ── generate teacher corpus ────────────────────────────────────────────────────
TEMP        = 0.6
GEN_LENGTH  = 300
NUM_SAMPLES = 500
EPOCHS      = 10
NGRAM_N     = 4

seeds = ["the ", "language ", "to be ", "a ", "in "]

generated_corpus = ""
for i in range(NUM_SAMPLES):
    generated_corpus += neural_generate(
        lstm_teacher, seeds[i % len(seeds)],
        char2id, id2char, id2idx, idx2id,
        vocab_size, char2cluster,
        length=GEN_LENGTH, temp=TEMP
    ) + " "

print(f"Teacher generated {len(generated_corpus)} chars of text.")

# ── build GPT-space vocab from generated corpus ────────────────────────────────
# Use word-level tokens (split on spaces) mapped into 0..VOCAB_SIZE-1
words  = generated_corpus.split()
unique = sorted(set(words))
# cap at VOCAB_SIZE
unique = unique[:VOCAB_SIZE]
w2i    = {w: i for i, w in enumerate(unique)}
i2w    = {i: w for w, i in w2i.items()}
gpt_vocab_size = len(unique)   # <= 50257; GPT net I/O stays 50257

# ── n-gram teacher over GPT token space ───────────────────────────────────────
ngram = NGramTeacher(n=NGRAM_N)
ngram.train(generated_corpus)   # char-level context, char-level counts

# soft target projected into word-token space:
# for each consecutive word pair, soft target = n-gram char distribution
# re-weighted by vortex math position scalar on the one-hot input

def vortex_weight(pos):
    """Scalar in (0,1] encoding vortex math position."""
    dr  = digital_root(pos)
    vp  = vortex_pos(dr)          # 0..5 or -1
    fg  = family_group(pos)       # 0,1,2
    w   = (vp + 1) / 6.0 if vp >= 0 else 0.5
    w  *= 1.0 - fg * 0.1          # family group nudge
    return max(w, 0.1)

def encode_token(token_id, pos):
    """One-hot over 50257, scaled by vortex weight at this position."""
    vec = np.zeros(VOCAB_SIZE)
    if token_id < VOCAB_SIZE:
        vec[token_id] = vortex_weight(pos)
    return vec

def ngram_soft_target_gpt(ctx_words, next_word_id):
    """
    Soft target over 50257: n-gram char distribution over the context,
    projected onto word ids by matching generated word chars.
    Falls back to hard one-hot if n-gram has no context.
    """
    ctx_str = ' '.join(ctx_words)
    char_dist = ngram.soft_target(ctx_str, sorted(char2id.keys()),
                                  {c: i for i, c in enumerate(sorted(char2id.keys()))},
                                  len(char2id))
    y = np.zeros(VOCAB_SIZE)
    if char_dist is not None:
        # distribute char probability mass to words whose first char matches
        for wid, word in i2w.items():
            if word and word[0] in char2id:
                cidx = sorted(char2id.keys()).index(word[0])
                y[wid] += char_dist[cidx]
        s = y.sum()
        if s > 0:
            y /= s
            return y
    # fallback: hard target
    y[next_word_id] = 1.0
    return y

# ── build dataset ──────────────────────────────────────────────────────────────
ds = SupervisedDataSet(VOCAB_SIZE, VOCAB_SIZE)
tokens = [w2i[w] for w in words if w in w2i]

CHUNK = 15
for start in range(0, len(tokens) - 1, CHUNK):
    seg = tokens[start:start + CHUNK + 1]
    for i in range(len(seg) - 1):
        pos      = start + i + 1
        x        = encode_token(seg[i], pos)
        ctx      = [i2w[t] for t in seg[max(0, i - NGRAM_N + 1):i + 1] if t in i2w]
        y        = ngram_soft_target_gpt(ctx, seg[i + 1])
        ds.addSample(x, y)

print(f"Distillation dataset: {len(ds)} samples  gpt_vocab={gpt_vocab_size}")

# ── load GPT student ───────────────────────────────────────────────────────────
from main_gpt3_3d import net
student = net

# ── distill ────────────────────────────────────────────────────────────────────
trainer = BackpropTrainer(student, ds, learningrate=0.0001, momentum=0.9)
for epoch in range(1, EPOCHS + 1):
    err = trainer.train()
    print(f"  epoch {epoch:3d}  error={err:.6f}")

# ── save ───────────────────────────────────────────────────────────────────────
NetworkWriter.writeToFile(student, "gpt_distilled.xml")
print("Saved gpt_distilled.xml")

# ── Test ───────────────────────────────────────────────────────────────────────
print("\n── GPT distilled responses ──────────────────────────────────────")

def gpt_generate(net, seed_words, length=20, temp=0.8):
    result = list(seed_words)
    for pos in range(len(seed_words), len(seed_words) + length):
        last = result[-1]
        tid  = w2i.get(last, 0)
        x    = encode_token(tid, pos)
        probs = np.array(net.activate(x), dtype=np.float64)
        # restrict to known vocab range
        probs[gpt_vocab_size:] = 0.0
        log_p = np.log(np.clip(probs, 1e-9, 1.0)) / temp
        log_p -= log_p.max()
        probs  = np.exp(log_p); probs /= probs.sum()
        next_id = np.random.choice(len(probs), p=probs)
        result.append(i2w.get(next_id, "<unk>"))
    return ' '.join(result)

for seed in [["the", "cat"], ["language", "is"], ["to", "be"]]:
    # filter seed to known vocab
    seed = [w for w in seed if w in w2i] or [words[0]]
    print(f"\n[seed='{' '.join(seed)}']")
    print(gpt_generate(student, seed, length=20, temp=0.8))

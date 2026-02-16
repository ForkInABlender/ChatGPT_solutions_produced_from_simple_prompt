# Dylan Kenneth Eliot

"""
Below is the code usable for bit-wise efficient layers.

Since not all layers need exact precision, carving down footprint was required.

What this allows is efficient use of memory. Given that this is to be adapted to the main__GPT_Model.py structure of GPT-3.5, it will further reduce overhead for the model.

"""


import numpy as np

def quantize(x, bits):
    if bits == 32:
        return x.astype(np.float32)
    if bits == 1:
        return np.sign(x)
    levels = 2 ** bits - 1
    max_val = np.max(np.abs(x)) + 1e-8
    scale = levels / max_val
    return np.round(x * scale) / scale

class MixedPrecisionModel:
    def __init__(self, vocab=32, d=16, hidden=32, seed=0):
        rng = np.random.default_rng(seed)
        self.embed = rng.standard_normal((vocab, d))
        self.fc1 = rng.standard_normal((d, hidden))
        self.fc2 = rng.standard_normal((hidden, d))
        self.out = rng.standard_normal((d, vocab))

        # mixed allocation
        self.embed_q = quantize(self.embed, 1)
        self.fc1_q   = quantize(self.fc1, 2)
        self.fc2_q   = quantize(self.fc2, 4)
        self.out_q   = quantize(self.out, 8)

    def forward(self, token_id):
        x = self.embed_q[token_id]
        x = np.tanh(x @ self.fc1_q)
        x = np.tanh(x @ self.fc2_q)
        return x @ self.out_q

class UniformModel:
    def __init__(self, bits, vocab=32, d=16, hidden=32, seed=0):
        rng = np.random.default_rng(seed)
        self.embed = rng.standard_normal((vocab, d))
        self.fc1 = rng.standard_normal((d, hidden))
        self.fc2 = rng.standard_normal((hidden, d))
        self.out = rng.standard_normal((d, vocab))

        self.embed_q = quantize(self.embed, bits)
        self.fc1_q   = quantize(self.fc1, bits)
        self.fc2_q   = quantize(self.fc2, bits)
        self.out_q   = quantize(self.out, bits)

    def forward(self, token_id):
        x = self.embed_q[token_id]
        x = np.tanh(x @ self.fc1_q)
        x = np.tanh(x @ self.fc2_q)
        return x @ self.out_q

# Run test
token = 5
seed = 42

mixed = MixedPrecisionModel(seed=seed)
uniform_bits = [1, 2, 4, 8, 16, 32]

results = {
    "with_mixed": {
        "pred": int(np.argmax(mixed.forward(token))),
        "mean": float(np.mean(mixed.forward(token))),
        "std": float(np.std(mixed.forward(token))),
        "top5": mixed.forward(token)[np.argsort(mixed.forward(token))[-5:]][::-1].round(3).tolist()
    }
}

for b in uniform_bits:
    u = UniformModel(bits=b, seed=seed)
    logits = u.forward(token)
    results[f"with_{b}bit"] = {
        "pred": int(np.argmax(logits)),
        "mean": float(np.mean(logits)),
        "std": float(np.std(logits)),
        "top5": logits[np.argsort(logits)[-5:]][::-1].round(3).tolist()
    }

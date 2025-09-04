# Dylan Kenneth Eliot

"""
This is how you GPT the wrong way, then pickle the object for layer reuse, minus the pickling of objects via modules such as 'dill' (https://pypi.org/project/dill/)

As you can see some inquiries might need adjustment of even the "weights" to come closer to matching quiry. The same might be also true for top-k, temp & top-p values.

This is the type of code I normally throw into the project-graveyard due to it's half-functional state.
"""

# Running a test of the offline GPT-3.5 behavioral emulator (trigram Markov + sampling)
# This cell defines the emulator (copied from the provided script) and runs several test prompts,
# saving outputs to ./emulator_test_results.txt and printing them here.

from collections import defaultdict, Counter
import random, math, re, sys

# --- Tokenize / detokenize utilities ---
def tokenize(text):
    text = text.strip()
    tokens = re.findall(r"\w+|[.,!?;:\'\"()\-\—]", text)
    return tokens

def detokenize(tokens):
    out = []
    prev = ""
    for t in tokens:
        if t.isalnum() or (len(t)==1 and t.isalpha()):
            if prev and (prev.isalnum() or prev in (')', '"', "'")):
                out.append(" ")
            out.append(t)
        else:
            out.append(t)
        prev = t
    s = "".join(out)
    s = re.sub(r"\s+([.,!?;:])", r"\1", s)
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r'(^|[.!?]\s+)([a-z])', lambda m: m.group(1) + m.group(2).upper(), s)
    return s

# --- NGramModel ---
class NGramModel:
    def __init__(self, order=3, smoothing=1e-3):
        self.order = order
        self.smoothing = smoothing
        self.counts = [defaultdict(Counter) for _ in range(order)]
        self.unigram = Counter()
        self.vocab = set()
        self.total_unigrams = 0

    def add_sentence(self, tokens):
        tokens = ["<s>"] + tokens + ["</s>"]
        for t in tokens:
            self.unigram[t] += 1
            self.vocab.add(t)
            self.total_unigrams += 1
        for n in range(2, self.order+1):
            for i in range(len(tokens)-n+1):
                prefix = tuple(tokens[i:i+n-1])
                nxt = tokens[i+n-1]
                self.counts[n-1][prefix][nxt] += 1

    def finalize(self):
        pass

    def candidate_probs(self, prefix):
        V = len(self.vocab) or 1
        smoothing = self.smoothing
        probs = {}
        prefix = tuple(prefix[-(self.order-1):]) if prefix else tuple()
        L = len(prefix)
        if L >= 2:
            weights = [0.1, 0.3, 0.6]
        elif L == 1:
            weights = [0.2, 0.8, 0.0]
        else:
            weights = [1.0, 0.0, 0.0]
        unigram_w, bigram_w, trigram_w = weights[0], weights[1], weights[2]
        uni_total = self.total_unigrams + smoothing * V
        uni_probs = {tok: (self.unigram.get(tok,0)+smoothing)/uni_total for tok in self.vocab}
        if L >= 1:
            big_prefix = (prefix[-1],)
            big_ctr = self.counts[1].get(big_prefix, Counter())
            big_total = sum(big_ctr.values()) + smoothing * V
            big_probs = {tok: (big_ctr.get(tok,0)+smoothing)/big_total for tok in self.vocab}
        else:
            big_probs = {tok: 0.0 for tok in self.vocab}
        if L >= 2:
            tri_prefix = (prefix[-2], prefix[-1])
            tri_ctr = self.counts[2].get(tri_prefix, Counter())
            tri_total = sum(tri_ctr.values()) + smoothing * V
            tri_probs = {tok: (tri_ctr.get(tok,0)+smoothing)/tri_total for tok in self.vocab}
        else:
            tri_probs = {tok: 0.0 for tok in self.vocab}
        for tok in self.vocab:
            probs[tok] = unigram_w * uni_probs.get(tok,0.0) + bigram_w * big_probs.get(tok,0.0) + trigram_w * tri_probs.get(tok,0.0)
        s = sum(probs.values()) or 1.0
        for k in list(probs.keys()):
            probs[k] = probs[k] / s
        return probs

# --- Sampling helpers ---
def apply_repetition_penalty(probs, recent_tokens, penalty=1.2):
    adjusted = {}
    recent_set = set(recent_tokens[-20:])
    for tok, p in probs.items():
        if tok in recent_set and p>0:
            adjusted[tok] = p / penalty
        else:
            adjusted[tok] = p
    s = sum(adjusted.values()) or 1.0
    for k in adjusted:
        adjusted[k] /= s
    return adjusted

def top_k_top_p_filter(probs, top_k=0, top_p=1.0):
    items = sorted(probs.items(), key=lambda x: x[1], reverse=True)
    if top_k > 0:
        items = items[:top_k]
    if top_p < 1.0:
        cum = 0.0
        kept = []
        for tok, p in items:
            kept.append((tok,p))
            cum += p
            if cum >= top_p:
                break
        items = kept
    s = sum(p for _, p in items) or 1.0
    return {tok: p/s for tok, p in items}

def sample_from_probs(probs, temperature=1.0):
    if temperature != 1.0 and temperature > 0:
        adjusted = {tok: math.pow(p, 1.0/temperature) for tok, p in probs.items()}
        s = sum(adjusted.values()) or 1.0
        adjusted = {tok: p/s for tok,p in adjusted.items()}
        probs = adjusted
    toks = list(probs.keys())
    weights = [probs[t] for t in toks]
    if not toks:
        return "</s>"
    return random.choices(toks, weights=weights, k=1)[0]

# --- Emulator class ---
class GPT35Emulator:
    def __init__(self, corpus_lines=None, order=3):
        self.model = NGramModel(order=order)
        if corpus_lines:
            for line in corpus_lines:
                toks = tokenize(line)
                if toks:
                    self.model.add_sentence(toks)
        else:
            builtin = [
                "A blockchain is a distributed ledger that records transactions across many computers.",
                "Water drips softly on tin roofs after a storm.",
                "Pros of remote work include flexibility and no commute; cons include loneliness and distraction.",
                "Turn off the water supply, disassemble the faucet, replace the washer, and reassemble.",
                "Captain Ahab's obsession with the white whale leads to tragedy in Moby Dick.",
                "Mix flour, milk, eggs and fry spoonfuls in butter until golden.",
                "Seasons are caused by Earth's axial tilt changing the angle of sunlight.",
                "Quantum entanglement links particles so that measuring one affects the other.",
                "Bean & Branch, Morning Roost, Steamy Pages, Copper Kettle, Dawn Brew.",
                "Dear Professor, I hope you are well. I am writing to request an extension due to unexpected circumstances."
            ]
            for line in builtin:
                toks = tokenize(line)
                self.model.add_sentence(toks)
        self.model.finalize()

    def generate(self, prompt, max_tokens=100, temperature=0.9, top_k=40, top_p=0.9, repetition_penalty=1.15, creativity="medium"):
        prompt_toks = tokenize(prompt)
        recent = ["<s>"] + prompt_toks[:]
        output = []
        for i in range(max_tokens):
            prefix = tuple(recent[-(self.model.order-1):]) if len(recent) > 0 else tuple()
            probs = self.model.candidate_probs(prefix)
            probs = apply_repetition_penalty(probs, recent, penalty=repetition_penalty)
            probs = top_k_top_p_filter(probs, top_k=top_k, top_p=top_p)
            tok = sample_from_probs(probs, temperature=temperature)
            if tok == "</s>":
                break
            recent.append(tok)
            output.append(tok)
            if len(output) > 8 and random.random() < 0.06:
                if output[-1] in (".", "?", "!"):
                    break
        text = detokenize(output)
        text = self.post_process_style(text, creativity)
        return text

    def post_process_style(self, text, creativity):
        if not text:
            return ""
        sentences = re.split(r'(?<=[.!?])\s+', text)
        if len(sentences) > 2 and len(max(sentences, key=len).split()) > 30:
            longest = max(sentences, key=len)
            clipped = " ".join(longest.split()[:25]) + "."
            sentences[sentences.index(longest)] = clipped
        text = " ".join(sentences).strip()
        if random.random() < 0.12:
            words = text.split()
            if len(words) > 4:
                i = random.randint(0, max(0, len(words)-4))
                phrase = " ".join(words[i:i+3])
                text = text + " " + phrase
        if re.search(r'\b(how to|explain|why|what)\b', text.lower()):
            if random.random() < 0.35:
                text = text.replace(".", ". It might help to try this approach.") if "." in text else text + " It might help to try this approach."
        if creativity == "low":
            words = text.split()
            if len(words) > 40:
                text = " ".join(words[:40]) + "..."
        elif creativity == "medium":
            words = text.split()
            if len(words) > 120:
                text = " ".join(words[:120]) + "..."
        return text

# --- Run tests ---
emulator = GPT35Emulator()

test_prompts = [
    "Explain blockchain simply.",
    "Write a haiku about rain.",
    "How to fix a leaky faucet?",
    "Summarize Moby Dick in one paragraph.",
    "Suggest five names for a coffee shop.",
    "Write a polite email to ask for a deadline extension."
]

results = []
for prompt in test_prompts:
    out_low = emulator.generate(prompt, max_tokens=80, temperature=0.6, top_k=30, top_p=0.85, creativity="low")
    out_med = emulator.generate(prompt, max_tokens=120, temperature=0.9, top_k=40, top_p=0.9, creativity="medium")
    out_high = emulator.generate(prompt, max_tokens=160, temperature=1.1, top_k=80, top_p=0.95, creativity="high")
    results.append((prompt, out_low, out_med, out_high))

# Save results to file and print summary
out_path = "./emulator_test_results.txt"
with open(out_path, "w", encoding="utf8") as fh:
    for prompt, low, med, high in results:
        fh.write(f"Prompt: {prompt}\n\nLOW creativity:\n{low}\n\nMEDIUM creativity:\n{med}\n\nHIGH creativity:\n{high}\n\n{'-'*60}\n\n")
print("Wrote test results to", out_path)
# Print results here for immediate feedback
for prompt, low, med, high in results:
    print("PROMPT:", prompt)
    print(" LOW:", low)
    print(" MED:", med)
    print(" HIGH:", high)
    print("-"*60)

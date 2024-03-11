# Dylan Kenneth Eliot & GPT-4-Plugins ( Alpha Edition )


"""

This reduces the complexity & allows adaptibility while abstracting the heavy lifting.

"""

import json
import os
import regex as re
from functools import lru_cache
import tiktoken

@lru_cache()
def bytes_to_unicode():
    """
    Returns a dictionary mapping all bytes to their corresponding unicode characters.
    """
    bs = (list(range(ord("!"), ord("~") + 1)) + list(range(ord("¡"), ord("¬") + 1)) +
          list(range(ord("®"), ord("ÿ") + 1)))
    cs = bs[:]
    n = 0
    for b in range(2**8):
        if b not in bs:
            bs.append(b)
            cs.append(2**8 + n)
            n += 1
    cs = [chr(n) for n in cs]
    return dict(zip(bs, cs))

def get_pairs(word):
    """
    Identifies all unique pairs of symbols (characters) in a word.
    """
    pairs = set()
    prev_char = word[0]
    for char in word[1:]:
        pairs.add((prev_char, char))
        prev_char = char
    return pairs

class Encoder:
    def __init__(self, model_name):
        # Initialize tiktoken encoding for the specified model
        self.encoding = tiktoken.get_encoding(model_name)
        self.byte_encoder = bytes_to_unicode()
        self.bpe_ranks = {}  # This would be populated based on the model's BPE merges

    def bpe(self, token):
        """
        Manually applies BPE to a token using the byte_encoder and bpe_ranks.
        """
        if token in self.bpe_ranks:
            return token  # Simplified; actual BPE logic is more complex

        word = tuple(self.byte_encoder[b] for b in token.encode('utf-8'))
        pairs = get_pairs(word)

        if not pairs:
            return token

        # Simplified BPE logic; actual implementation would use bpe_ranks to merge pairs
        while True:
            bigram = min(pairs, key=lambda pair: self.bpe_ranks.get(pair, float('inf')))
            if bigram not in self.bpe_ranks:
                break
            # Actual merging logic would go here
            pairs.remove(bigram)

        # This is a placeholder to indicate where the BPE merge logic would be applied
        return token

    def encode(self, text):
        """
        Encodes text using a combination of manual BPE and tiktoken encoding.
        """
        encoded_tokens = []
        for token in text.split():  # Simplified tokenization
            bpe_token = self.bpe(token)
            encoded_tokens.extend(self.encoding.encode(bpe_token))
        return encoded_tokens

    def decode(self, tokens):
        """
        Decodes tokens to text using tiktoken's decoding.
        """
        return self.encoding.decode(tokens)

def get_encoder(model_name):
    # Initialize the Encoder with the specified model_name
    # This function acts as a factory to create Encoder instances
    return Encoder(model_name)

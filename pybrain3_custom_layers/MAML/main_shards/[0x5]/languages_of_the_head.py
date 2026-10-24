
# Dylan Kenneth Eliot
"""
Context-dependent multilingual compact encoder and decoder.

Encoding routes:
    1. Known English words -> Hindi vocabulary.
    2. Other known English words -> Chinese vocabulary.
    3. Unmapped words -> opaque glyph sequences.

Features:
    - Python 3.10+, standard library only.
    - Offline and deterministic within an encoder session.
    - Preserves whitespace, punctuation, and line breaks.
    - Supports multiword translations.
    - Preserves original spelling and capitalization.
    - Includes a reversible decoder.
    - No NLTK, TextBlob, or network downloads.

Note:
    This is a reversible encoding as encryption.
    The glyph vocabulary must remain available for decoding.
    If you manually change the order id each vocabulary defined in the file uses, remember that 0 has to exist.
     and negative numbers do not aide in translation.
"""

import re
from typing import Dict, List, Tuple


# ------------------------------------------------------------
# Glyph inventory
# ------------------------------------------------------------

SYMBOLS: Tuple[str, ...] = tuple("""
𓂀 𓆣 𓂸 𓅓 𓃰 𓀠 𓃱 𓂻 𓆏 𓀿 𓁿 𓃬 𓂃
𓃒 𓀒 𓃠 𓆙 ⚆ ╯ ┻ ━ ╥ 益 ಥ ⌑ ⌒ °
□ ´ ； ˶ ∩ 彡 ⁀ ⊙ ☉ º Δ ༎ ຶ ⇀
‸ T ⌓ ╬ 𓄿 𓅂 𓆎 ⚘ ⍨ ᚠ ᚢ ᚦ ᚨ ᚱ
ᚲ ᚷ ᚹ ᚺ ᚾ ᛁ ᛃ ᛈ ᛇ ᛉ ᛋ ᛏ ᛒ ᛖ
ᛗ ᛚ ᛜ ᛞ
""".split())

BASE = len(SYMBOLS)

if BASE < 2:
    raise RuntimeError("At least two unique symbols are required.")

if len(set(SYMBOLS)) != BASE:
    raise RuntimeError("The glyph inventory contains duplicate symbols.")


# ------------------------------------------------------------
# Translation dictionaries
# ------------------------------------------------------------

HINDI_TRANSLATION: Dict[str, str] = {
    "the": "यह",
    "of": "का",
    "and": "और",
    "in": "में",
    "to": "को",
    "a": "अ",
    "is": "है",
    "for": "के लिए",
    "that": "वह",
    "on": "पर",
    "with": "साथ",
    "as": "जैसा",
    "by": "द्वारा",
    "an": "अन",
    "be": "होना",
    "are": "हैं",
    "at": "पर",
    "from": "से",
    "this": "यह",
    "or": "या",
    "but": "लेकिन",
    "not": "नहीं",
}

CHINESE_TRANSLATION: Dict[str, str] = {
    "hello": "你好",
    "world": "世界",
    "memory": "记忆",
    "information": "信息",
    "language": "语言",
    "knowledge": "知识",
    "time": "时间",
    "person": "人",
    "people": "人们",
    "meaning": "意义",
    "context": "上下文",
    "word": "词",
    "words": "词语",
    "write": "写",
    "text": "文本",
    "long": "长",
    "short": "短",
    "unknown": "未知",
    "true": "真",
    "false": "假",
}


# ------------------------------------------------------------
# Tokenization
# ------------------------------------------------------------

# Whitespace, words, numbers, and punctuation are separate tokens.
# Concatenating the tokens reconstructs the original text.
TOKEN_PATTERN = re.compile(
    r"\s+"
    r"|[^\W\d_]+(?:['’][^\W\d_]+)*"
    r"|\d+(?:[.,]\d+)*"
    r"|[^\w\s]"
    r"|\w+",
    re.UNICODE,
)

WORD_PATTERN = re.compile(
    r"^[^\W\d_]+(?:['’][^\W\d_]+)*$",
    re.UNICODE,
)


# ------------------------------------------------------------
# Encoder and decoder
# ------------------------------------------------------------

import re
import json
from pathlib import Path



class CompactEncoder:
    """Reversible encoder without open, close, or separator markers."""

    def __init__(self):
        self.word_to_id = {}
        self.id_to_word = {}

        self.symbols = tuple(SYMBOLS)
        self.base = len(self.symbols)

        if self.base < 2:
            raise ValueError("At least two symbols are required.")

        if len(set(self.symbols)) != self.base:
            raise ValueError("SYMBOLS contains duplicate glyphs.")

        self.glyph_to_index = {
            glyph: index
            for index, glyph in enumerate(self.symbols)
        }

        # Each word ID uses the same number of glyphs.
        self.width = 1

    def _update_width(self):
        """Ensure the fixed-width representation fits all word IDs."""
        capacity = self.base

        while capacity < max(1, len(self.word_to_id)):
            self.width += 1
            capacity *= self.base

    def encode_id(self, number):
        """Encode an integer using exactly self.width glyphs."""
        if number < 0:
            raise ValueError("number must be nonnegative")

        digits = []

        for _ in range(self.width):
            number, remainder = divmod(number, self.base)
            digits.append(self.symbols[remainder])

        if number:
            raise ValueError("Glyph ID exceeds current encoding width.")

        return "".join(reversed(digits))

    def decode_id(self, glyph_sequence):
        """Decode a fixed-width glyph sequence."""
        number = 0

        for glyph in glyph_sequence:
            if glyph not in self.glyph_to_index:
                raise ValueError(f"Unknown glyph: {glyph!r}")

            number = (
                number * self.base
                + self.glyph_to_index[glyph]
            )

        return number

    def translate_word(self, word):
        """Use the original dictionaries, then glyph-encode unknown words."""
        key = word.casefold()

        if key in HINDI_TRANSLATION:
            return HINDI_TRANSLATION[key]

        if key in CHINESE_TRANSLATION:
            return CHINESE_TRANSLATION[key]

        return "".join(
            self.symbols[ord(char.lower()) - ord("a")]
            for char in word
        )

    def encode(self, text):
        """
        Encode words directly as fixed-width glyph IDs.

        Whitespace and punctuation are retained literally.
        """
        output = []
        position = 0

        for match in re.finditer(r"[A-Za-z]+", text):
            output.append(text[position:match.start()])
            word = match.group(0)

            if word not in self.word_to_id:
                index = len(self.word_to_id)
                self.word_to_id[word] = index
                self.id_to_word[index] = word

            self._update_width()
            output.append(self.encode_id(self.word_to_id[word]))
            position = match.end()

        output.append(text[position:])
        return "".join(output)

    def decode(self, encoded_text):
        """
        Decode fixed-width glyph sequences using the saved vocabulary.

        Non-glyph characters are preserved as literal text.
        """
        output = []
        glyph_run = []

        def flush():
            if not glyph_run:
                return

            sequence = "".join(glyph_run)
            glyph_run.clear()

            if len(sequence) % self.width:
                raise ValueError("Incomplete glyph ID sequence.")

            for offset in range(0, len(sequence), self.width):
                chunk = sequence[offset:offset + self.width]
                word_id = self.decode_id(chunk)

                if word_id not in self.id_to_word:
                    raise ValueError(f"Unknown word ID: {word_id}")

                output.append(self.id_to_word[word_id])

        for char in encoded_text:
            if char in self.glyph_to_index:
                glyph_run.append(char)
            else:
                flush()
                output.append(char)

        flush()
        return "".join(output)


        def restore(match):
            glyph_id = self.decode_id(match.group(1))

            if glyph_id not in self.id_to_word:
                raise ValueError(
                    f"Unknown word ID: {glyph_id}"
                )

            return self.id_to_word[glyph_id]

        decoded = pattern.sub(restore, encoded_text)

        if self.OPEN in decoded or self.CLOSE in decoded:
            raise ValueError("Malformed encoded text.")

        return decoded

    def save_vocabulary(self, filename):
        """Persist the session vocabulary for later decoding."""
        with open(filename, "w", encoding="utf-8") as file:
            json.dump(
                {
                    "word_to_id": self.word_to_id,
                    "id_to_word": {
                        str(key): value
                        for key, value in self.id_to_word.items()
                    },
                },
                file,
                ensure_ascii=False,
                indent=2,
            )

    def load_vocabulary(self, filename):
        """Restore a previously saved session vocabulary."""
        with open(filename, "r", encoding="utf-8") as file:
            data = json.load(file)

        self.word_to_id = {
            key: int(value)
            for key, value in data["word_to_id"].items()
        }

        self.id_to_word = {
            int(key): value
            for key, value in data["id_to_word"].items()
        }



# ------------------------------------------------------------
# Convenience function
# ------------------------------------------------------------

def translate_full_sentence(sentence: str) -> str:
    """
    Encode one sentence using a new session.

    Keep a CompactEncoder instance if you also need to decode
    the result later.
    """
    return CompactEncoder().encode(sentence)


# ------------------------------------------------------------
# Example
# ------------------------------------------------------------

encoder = CompactEncoder()

#original = (
#    "The memory stores information. Hello world!\n"
#    "abcdefghijklmnopqrstuvwxyz\n"
#    "zzzzzz aaaa abcd qzxv"
#)

#encoded = encoder.encode(original)
#encoder.save_vocabulary("language_vocabulary.json")

encoded = (
 "𓂀 𓆣 𓂸 𓅓. 𓃰 𓀠!\n"
"𓃱\n"
"𓂻 𓆏 𓀿 𓁿"
)

encoder.load_vocabulary("language_vocabulary.json")

decoded = encoder.decode(encoded)

#print("Original:", original)
print("Encoded:", encoded)
print("Decoded:", decoded)

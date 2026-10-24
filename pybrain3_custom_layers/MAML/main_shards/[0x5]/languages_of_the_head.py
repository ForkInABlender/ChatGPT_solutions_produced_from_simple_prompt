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
"""

import re
import json
from typing import Dict, Tuple


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
# Encoder and decoder
# ------------------------------------------------------------

class CompactEncoder:
    """Context-dependent multilingual compact encoder and decoder."""

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

        # Reverse dictionaries for decoding phrase lookups
        self.hindi_to_english = {v: k for k, v in HINDI_TRANSLATION.items()}
        self.chinese_to_english = {v: k for k, v in CHINESE_TRANSLATION.items()}

        # Each word ID uses a fixed number of glyphs.
        self.width = 1

    def _update_width(self):
        """Ensure the fixed-width representation fits all word IDs."""
        capacity = self.base

        while capacity < max(1, len(self.word_to_id)):
            self.width += 1
            capacity *= self.base

    def encode_id(self, number: int) -> str:
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

    def decode_id(self, glyph_sequence: str) -> int:
        """Decode a fixed-width glyph sequence into an integer ID."""
        number = 0

        for glyph in glyph_sequence:
            if glyph not in self.glyph_to_index:
                raise ValueError(f"Unknown glyph: {glyph!r}")

            number = (
                number * self.base
                + self.glyph_to_index[glyph]
            )

        return number

    def encode(self, text: str) -> str:
        """
        Encode text using:
            1. Hindi dictionary for matching words.
            2. Chinese dictionary for matching words.
            3. Fixed-width glyph sequences for unmapped words.
        Whitespace and punctuation are retained literally.
        """
        output = []
        position = 0

        for match in re.finditer(r"[A-Za-z]+", text):
            output.append(text[position:match.start()])
            word = match.group(0)
            key = word.casefold()

            if key in HINDI_TRANSLATION:
                output.append(HINDI_TRANSLATION[key])
            elif key in CHINESE_TRANSLATION:
                output.append(CHINESE_TRANSLATION[key])
            else:
                if word not in self.word_to_id:
                    index = len(self.word_to_id)
                    self.word_to_id[word] = index
                    self.id_to_word[index] = word

                self._update_width()
                output.append(self.encode_id(self.word_to_id[word]))

            position = match.end()

        output.append(text[position:])
        return "".join(output)

    def decode(self, encoded_text: str) -> str:
        """
        Decode encoded text back to original English words using reverse Hindi,
        Chinese dictionaries, and glyph ID sequences.
        """
        # Sort phrases by length descending to support greedy multi-word matching
        hindi_phrases = sorted(self.hindi_to_english.keys(), key=len, reverse=True)
        chinese_phrases = sorted(self.chinese_to_english.keys(), key=len, reverse=True)

        output = []
        i = 0
        n = len(encoded_text)

        while i < n:
            # 1. Check for glyph sequences (unmapped words)
            if encoded_text[i] in self.glyph_to_index:
                glyph_run = []
                while i < n and encoded_text[i] in self.glyph_to_index:
                    glyph_run.append(encoded_text[i])
                    i += 1

                sequence = "".join(glyph_run)
                if len(sequence) % self.width:
                    raise ValueError("Incomplete glyph ID sequence.")

                for offset in range(0, len(sequence), self.width):
                    chunk = sequence[offset:offset + self.width]
                    word_id = self.decode_id(chunk)

                    if word_id not in self.id_to_word:
                        raise ValueError(f"Unknown word ID: {word_id}")

                    output.append(self.id_to_word[word_id])
                continue

            # 2. Check for Hindi phrase match
            matched = False
            for phrase in hindi_phrases:
                if encoded_text.startswith(phrase, i):
                    output.append(self.hindi_to_english[phrase])
                    i += len(phrase)
                    matched = True
                    break
            if matched:
                continue

            # 3. Check for Chinese phrase match
            for phrase in chinese_phrases:
                if encoded_text.startswith(phrase, i):
                    output.append(self.chinese_to_english[phrase])
                    i += len(phrase)
                    matched = True
                    break
            if matched:
                continue

            # 4. Preserve literal characters (whitespace, punctuation, etc.)
            output.append(encoded_text[i])
            i += 1

        return "".join(output)

    def save_vocabulary(self, filename: str):
        """Persist the session vocabulary for later decoding."""
        with open(filename, "w", encoding="utf-8") as file:
            json.dump(
                {
                    "width": self.width,
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

    def load_vocabulary(self, filename: str):
        """Restore a previously saved session vocabulary."""
        with open(filename, "r", encoding="utf-8") as file:
            data = json.load(file)

        self.width = data.get("width", 1)
        self.word_to_id = {
            key: int(value)
            for key, value in data["word_to_id"].items()
        }

        self.id_to_word = {
            int(key): value
            for key, value in data["id_to_word"].items()
        }



# ------------------------------------------------------------
# Example
# ------------------------------------------------------------

encoder = CompactEncoder()

original = (
"The memory stores information. Hello world!\n"
"UnmappedWord example."
)
encoded = encoder.encode(original)
encoder.save_vocabulary("language_vocabulary.json")

# Decoder instance (loaded with vocabulary for unmapped words)
decoder = CompactEncoder()
decoder.load_vocabulary("language_vocabulary.json")
decoded = decoder.decode(encoded)

print("Original:", original)
print("\nEncoded:", encoded)
print("\nDecoded:", decoded)

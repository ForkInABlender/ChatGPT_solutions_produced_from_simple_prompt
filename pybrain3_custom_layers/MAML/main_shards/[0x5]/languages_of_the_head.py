# Dylan Kenneth Eliot

"""
This language is a language of languages.

It contains:
 - heiroglyphs
 - hindi
 - chinese simplified

in an order knowing any one language will not give the proper translation of any one part.
This language is a combination of the 3 languages. So, you'll have to use this to decode.

This scripts' purpose is for informational compactness. Meaning all translation is done to condense information based on how we'd store it later for long term memory but accessible via short term memory recall to long
 term memory store.

This informational compactness does encode and intentionally I have left off the decoder. This was so as to keep the information known to those that know the writer of the text and the context both know of, for which
 it was encoded. Either they know it inherently, or they do not and have to rely on a really good guess. If you don't have the context they share with the initial writer, it remains unknown to the ones translating.

"""

import nltk
from nltk.tokenize import word_tokenize

# List of 72 symbols to be used for translation
symbols = [
    "𓂀", "𓆣", "𓂸", "𓅓", "𓃰", "𓀠", "𓃱", "𓂻", "𓆏", "𓀿", "𓁿", "𓃬", "𓂃",
    "𓃒", "𓀒", "𓃠", "𓆙", "⚆", "╯", "┻", "━", "╥", "益", "ಥ", "⌑", "⌒", "°",
    "□", "´", "；", "˶", "∩", "彡", "⁀", "⊙", "☉", "º", "Δ", "༎", "ຶ", "⇀", 
    "‸", "T", "⌓", "╬", "𓄿", "𓅂", "𓆎", "⚘", "⍨", "ᚠ", "ᚢ", "ᚦ", "ᚨ", "ᚱ",
    "ᚲ", "ᚷ", "ᚹ", "ᚺ", "ᚾ", "ᛁ", "ᛃ", "ᛈ", "ᛇ", "ᛉ", "ᛋ", "ᛏ", "ᛒ", "ᛖ",
    "ᛗ", "ᛚ", "ᛜ", "ᛞ"
]

# Dictionary for common words translated to Hindi
hindi_translation = {
    'the': 'यह',
    'of': 'का',
    'and': 'और',
    'in': 'में',
    'to': 'को',
    'a': 'एक',
    'is': 'है',
    'for': 'के लिए',
    'that': 'वह',
    'on': 'पर',
    'with': 'साथ',
    'as': 'जैसा',
    'by': 'द्वारा',
    'an': 'एक',
    'be': 'होना',
    'are': 'हैं',
    'at': 'पर',
    'from': 'से',
    'this': 'यह',
    'or': 'या',
    'but': 'लेकिन',
    'not': 'नहीं'
}

def translate_word(word, index):
    """Translates a word based on the Hindi dictionary or assigns a symbol if the word is not found."""
    # If the word is in the Hindi dictionary, translate it
    if word.lower() in hindi_translation:
        return hindi_translation[word.lower()]
    # Otherwise, use the symbols for translation, cycling through the list
    elif word.isalnum():
        return symbols[index % len(symbols)]
    # Leave punctuation intact
    return word

def translate_full_sentence(sentence):
    """Translates each word of a sentence to either a Hindi translation or a unique symbol."""
    # Tokenize the sentence using NLTK
    words = word_tokenize(sentence)

    # Translate each word individually using either Hindi dictionary or symbols
    translated_words = [translate_word(word, index) for index, word in enumerate(words)]

    # Reconstruct the translated sentence
    translated_sentence = ' '.join(translated_words)
    return translated_sentence

# Example usage

sentence = "No more guns."
translated_sentence = translate_full_sentence(sentence)

print("Original Sentence:", sentence)
print("Translated Sentence:", translated_sentence)

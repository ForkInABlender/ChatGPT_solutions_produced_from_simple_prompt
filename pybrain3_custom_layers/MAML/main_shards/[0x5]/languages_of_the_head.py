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


"""

# Required Libraries
import re
import nltk
from nltk.corpus import wordnet, stopwords
from collections import defaultdict

# Given Text Data (same as before)
data = '''
𓂀 (Eye of Horus) - Protection, spiritual insight, or vision.
Usage: Denote awareness or seeking divine wisdom, especially when invoking a metaphorical "seeing the truth."

𓆣 (Lotus Flower) - Growth, purity, beauty.
Usage: Convey themes of rebirth, enlightenment, or pure intentions.

𓂸 (Depiction of a Man Kneeling) - Humility or respect.
Usage: Illustrate moments of surrender, bowing down to wisdom, or when asking for guidance.

𓅓 (Ibis Bird) - Knowledge, Thoth's wisdom.
Usage: Use in educational contexts or when seeking or providing deep insight.

𓃰 (Ram's Head) - Strength, virility.
Usage: Represent resilience or power, particularly in challenging situations.

𓀠 (Standing Man) - Courage, presence.
Usage: Demonstrate the need to face adversity or stand one’s ground.

𓃱 (Hippopotamus) - Protection, maternal instinct.
Usage: Invoke in contexts of defending loved ones or safeguarding something valuable.

𓂻 (Foot or Leg in Motion) - Progress, movement.
Usage: Illustrate journeys, either physical or mental, or the process of moving forward.

𓆏 (Frog) - Fertility, transformation.
Usage: Denote changes or transformations, particularly with themes of nature and growth.

𓀿 (Pointing Hand) - Guidance, focus.
Usage: Direct attention or give guidance—especially when trying to make a point clear.

𓁿 (Lion) - Strength, royalty.
Usage: Address leadership or courage in facing challenges.

𓃬 (Crocodile) - Power, hidden danger.
Usage: Discuss potential threats or challenges lying beneath the surface.

𓂃 (Outstretched Arm) - Helping, action.
Usage: Emphasize aid or extending help to others.

𓃒 (Cat) - Protection, grace.
Usage: Emphasize independence or watchfulness.

𓀒 (Sitting Man) - Patience, contemplation.
Usage: Denote the importance of taking time for reflection or waiting.

𓃠 (Dog) - Loyalty, companionship.
Usage: Denote the value of friendship or trust.

𓆙 (Snake) - Renewal, healing.
Usage: Employ in contexts of transformation, particularly related to self-renewal or overcoming obstacles.

⚆ (Circle with Eyes) - Curiosity, surprise.
Usage: Express being intrigued or bewildered by something.

╯ (Left Parenthesis for Movement) - Action, motion.
Usage: Indicate active gestures, particularly related to throwing or shifting objects.

┻ (Falling Over Table) - Anger, frustration.
Usage: Use in chaotic situations or when expressing frustration.

━ (Straight Line for Barriers) - Boundaries, limits.
Usage: Define limitations or set distinct borders.

╥ (Fountain of Tears) - Sadness, release.
Usage: Use during emotional outbursts or shedding of tears.

益 (Japanese 'Benefit') - Growth, positivity.
Usage: Show progression or the benefit of a difficult situation.

ಥ (Tearful Face) - Emotional impact.
Usage: Illustrate hurt feelings or vulnerability.

⌑ (Symbol for Geometry) - Calculation, thoughtfulness.
Usage: Emphasize the need for precision or thought.

⌒ (Curved Line) - Fluidity, flexibility.
Usage: Discuss adapting or bending with ease in different contexts.

° (Small Circle) - Completeness.
Usage: Denote the smallness of an element in comparison or show a complete cycle.

□ (Square) - Structure, stability.
Usage: Symbolize rigidity or organization.

´ (Accent Mark) - Emphasis.
Usage: Stress a particular sound or detail.

； (Semicolon) - Continuation.
Usage: Bridge two closely related ideas.

˶ (Quotation-like Mark) - Emphasis on speech or action.
Usage: Frame quoted or emphasized words.

∩ (Intersection) - Connection or overlap.
Usage: Denote a connection or shared point between two things.

彡 (Three Small Streaks) - Movement or energy.
Usage: Indicate swift movement, often in a humorous or dynamic context.

⁀ (Curved Line for Smile) - Positivity, hope.
Usage: Indicate a hopeful or light-hearted approach.

⊙ (Circle with Dot) - Clarity or focus.
Usage: Convey an enlightened or clearly understood idea.

☉ (Sun Symbol) - Light, guidance.
Usage: Invoke clarity, warmth, or the concept of illuminating the truth.

º (Degree Symbol) - Temperature, measurement.
Usage: Show measurement or specific value.

Δ (Delta, Triangle) - Change, transformation.
Usage: Signify shifts or transformations in thought or action.

༎ (Teardrop-like Symbol) - Sorrow, heaviness.
Usage: Illustrate situations involving deep emotion or sadness.

ຶ (Small Mark/Accent) - Subtlety.
Usage: Show small but significant shifts.

⇀ (Arrow Right) - Moving Forward.
Usage: Denote direction or motivation.

‸ (Caret-like Symbol) - Highlighting.
Usage: Indicate importance, typically underlining what comes next.

T (Letter T) - Signify support.
Usage: Indicate a crossroads or decision point.

⌓ (Rounded Square) - Confusion, puzzlement.
Usage: Express a mixed or complicated feeling.

╬ (Cross with Multiple Lines) - Complex interrelations.
Usage: Denote the intersection of multiple ideas or paths.

⚘ (Flower) - Peace, beauty.
Usage: Denote growth, gifting, or peaceful endeavors.

⍨ (Double Arch) - Unease, playfulness.
Usage: Convey a quirky or uncertain tone.
'''

# Function to parse the data and build the phrase map (same as before)
def build_phrase_map(data):
    entries = data.strip().split('\n\n')
    phrase_map = {}
    for entry in entries:
        lines = entry.strip().split('\n')
        if len(lines) >= 2:
            # Extract symbol and meaning
            symbol_line = lines[0]
            usage_line = lines[1]

            # Parse symbol and meaning
            try:
                symbol_part, meaning_part = symbol_line.split(') - ')
                symbol = symbol_part.strip().split(' (')[0]
                meaning = meaning_part.strip()
            except ValueError:
                continue  # Skip malformed lines

            # Parse usage
            usage = usage_line.replace('Usage:', '').strip()

            # Add to phrase map
            phrase_map[symbol] = {'meaning': meaning, 'usage': usage}
    return phrase_map

# Function to build a word-to-symbol mapping with priorities (same as before)
def build_word_symbol_map(phrase_map):
    word_symbol_map = defaultdict(lambda: {})
    stop_words = set(stopwords.words('english'))
    for symbol, info in phrase_map.items():
        # Process meaning with priority 3
        meaning_text = info['meaning']
        meaning_words = re.findall(r'\b\w+\b', meaning_text.lower())
        for word in meaning_words:
            if word in stop_words:
                continue
            word_symbol_map[word][symbol] = 3  # Highest priority

        # Process usage with priority 2
        usage_text = info['usage']
        usage_words = re.findall(r'\b\w+\b', usage_text.lower())
        for word in usage_words:
            if word in stop_words:
                continue
            if symbol not in word_symbol_map[word]:
                word_symbol_map[word][symbol] = 2

        # Process synonyms with priority 1
        combined_text = meaning_text + ' ' + usage_text
        words = re.findall(r'\b\w+\b', combined_text.lower())
        for word in words:
            if word in stop_words:
                continue
            synsets = wordnet.synsets(word)
            for syn in synsets:
                for lemma in syn.lemmas():
                    lemma_name = lemma.name().replace('_', ' ').lower()
                    if lemma_name in stop_words:
                        continue
                    if symbol not in word_symbol_map[lemma_name]:
                        word_symbol_map[lemma_name][symbol] = 1
    return word_symbol_map

# Manual mapping of English stopwords to Hindi (same as before)
hindi_stopwords = {
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
    'not': 'नहीं',
    'requires': 'आवश्यक है',
    'here':"यहाँ",
    'there':"वहाँ",
    'their':"उनका",
    "they": "वे।",
    'perhaps':"शायद",
    # Add more stopwords as needed
}

# Function to translate words to Chinese using NLTK's WordNet
def translate_to_chinese(word):
    synsets = wordnet.synsets(word, lang='eng')
    for syn in synsets:
        # Try to get Chinese lemmas
        chinese_lemmas = syn.lemma_names('cmn')
        if chinese_lemmas:
            # Return the first Chinese lemma found
            return chinese_lemmas[0]
    # If no Chinese translation found, return the original word
    return word

# Function to translate text to the new symbol set, translate stopwords to Hindi, and remaining words to Chinese
def translate_to_symbols(text, word_symbol_map, hindi_stopwords):
    # Tokenize the input text
    words = re.findall(r'\b\w+\b|\S', text)
    translated = []
    for word in words:
        key = word.lower()
        if key in word_symbol_map:
            # Select the symbol with the highest priority
            symbols_with_priorities = word_symbol_map[key]
            sorted_symbols = sorted(symbols_with_priorities.items(), key=lambda x: x[1], reverse=True)
            symbol = sorted_symbols[0][0]
            translated.append(symbol)
        elif key in hindi_stopwords:
            # Replace stopword with its Hindi translation
            translated.append(hindi_stopwords[key])
        else:
            # Translate the word to Chinese
            chinese_word = translate_to_chinese(key)
            translated.append(chinese_word)
    return ''.join(translated)

# Build the phrase map from the data
phrase_map = build_phrase_map(data)

# Build the word-to-symbol mapping with priorities
word_symbol_map = build_word_symbol_map(phrase_map)

# Example usage of the translate_to_symbols function
input_text = "input text goes there. The rest will remain connected. Just fill in the blanks. Their is likely a full translation even for the couch. Perhaps they're as they are."
translated_text = translate_to_symbols(input_text, word_symbol_map, hindi_stopwords)
print("Original Text:")
print(input_text)
print("\nTranslated Text:")
print(translated_text)

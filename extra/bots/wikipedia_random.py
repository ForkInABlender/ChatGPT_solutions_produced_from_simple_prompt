# Dylan Kenneth Eliot

"""
True random is comes from controlled orders of chaotic rambunction following normal rules of logic while displaying absolute disobedience to what should be the
 rule.
"""

import wikipedia
import random
import re

# Word banks
determiners = ["the", "a", "every", "some", "each"]
nouns = ["cat", "robot", "ocean", "banana", "flame", "vortex", "pickle", "shadow", "toaster", "quantum"]
verbs = ["eats", "jumps", "vaporizes", "sings", "melts", "copies", "chases", "teleports", "encrypts", "shuffles"]
adjectives = ["blue", "fiery", "invisible", "slippery", "cosmic", "dizzy", "feral", "synthetic", "quiet", "jagged"]
prepositions = ["under", "over", "beside", "into", "through", "across", "beyond", "within", "near", "against"]
conjunctions = ["and", "but", "or", "so", "yet"]

# Function to get a random Wikipedia article
def get_random_wikipedia_content():
    try:
        # Fetch a random Wikipedia page
        page_title = wikipedia.random()
        
        # Search if the page is a disambiguation page and handle it
        search_results = wikipedia.search(page_title)
        
        if len(search_results) > 1:
            # If disambiguation, take the first valid result that isn't a disambiguation page
            for result in search_results:
                try:
                    # Try to fetch the page for each search result
                    page = wikipedia.page(result)
                    if page.content:
                        print(f"Using page: {result}")
                        return page.content
                except wikipedia.exceptions.DisambiguationError:
                    continue
            print(f"Disambiguation page found for: {page_title}")
            return ""
        else:
            # If it's not disambiguated, fetch the content normally
            content = wikipedia.page(page_title).content
            return content
    except wikipedia.exceptions.DisambiguationError as e:
        print(f"DisambiguationError: {e}")
        return ""
    except Exception as e:
        print(f"Error fetching Wikipedia content: {e}")
        return ""

# Clean content to extract words for filling template
def extract_words_from_content(content):
    # Remove non-alphanumeric characters and split into words
    words = re.findall(r'\b\w+\b', content.lower())
    return words

# Function to fill the template with random Wikipedia content
def generate_sentence_from_wikipedia(pattern: str, count=1):
    pattern_tags = pattern.strip().split()
    results = []

    for _ in range(count):
        sentence = []
        # Get random Wikipedia content
        wikipedia_content = get_random_wikipedia_content()
        if not wikipedia_content:
            # If no content fetched, fall back to random words
            print("Falling back to random words due to content issue.")
            wikipedia_content = " ".join(random.choices(nouns + verbs + adjectives, k=50))

        words = extract_words_from_content(wikipedia_content)

        for tag in pattern_tags:
            if tag == "[Determiner]":
                sentence.append(random.choice(determiners))
            elif tag == "[Noun]":
                sentence.append(random.choice(nouns) if not words else random.choice(words))
            elif tag == "[Verb]":
                sentence.append(random.choice(verbs) if not words else random.choice(words))
            elif tag == "[Adjective]":
                sentence.append(random.choice(adjectives) if not words else random.choice(words))
            elif tag == "[Preposition]":
                sentence.append(random.choice(prepositions))
            elif tag == "[Conjunction]":
                sentence.append(random.choice(conjunctions))
            else:
                sentence.append("[UNKNOWN]")

        results.append(" ".join(sentence))

    return results

# Example pattern usage
pattern = "[Determiner] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Adjective] [Noun] [Preposition] [Adjective] [Noun] [Conjunction] [Verb] [Adjective] [Noun] [Preposition] [Adjective] [Noun] [Preposition] [Determiner] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Conjunction] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective] [Noun] [Verb] [Determiner] [Adjective] [Noun] [Preposition] [Determiner] [Adjective]"
for sentence in generate_sentence_from_wikipedia(pattern, count=5):
    print(sentence)

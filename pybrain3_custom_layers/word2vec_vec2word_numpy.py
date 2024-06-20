# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is useful with pybrain3 neural networks that have the correct parameters and interfacing. 

This removes the need for tiktoken api offline. 
(https://github.com/openai/tiktoken & https://stackoverflow.com/questions/76106366/how-to-use-tiktoken-in-offline-mode-computer for point of reference) 

"""

import numpy as np

def sentence_to_ascii_sublists(sentence):
    words = sentence.split(' ')
    ascii_sublists = [[ord(char) for char in word] for word in words]
    return ascii_sublists

def ascii_sublists_to_sentence(ascii_sublists):
    words = [''.join(chr(value) for value in sublist) for sublist in ascii_sublists]
    return ' '.join(words)

sentence = "Using NumPy to parse words in a sentence is interesting."
ascii_sublists = sentence_to_ascii_sublists(sentence)
print("Original Sentence:", sentence)
print("ASCII Sublists:", ascii_sublists)
print("\nReconstructed Sentence from ASCII Sublists:", ascii_sublists_to_sentence(ascii_sublists))

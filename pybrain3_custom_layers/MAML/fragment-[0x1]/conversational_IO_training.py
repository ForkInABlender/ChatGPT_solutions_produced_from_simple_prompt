# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is the basic template that will be refactored for more in-depth AI development. 

Because it expects inputs and outputs to be 1-d arrays, all input and output to the model should be a 1-d arrays.

Later on, this will be useful with the dim3_neuronlayer class for improved informational processing.

https://github.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/blob/2024_06/pybrain3_custom_layers/dim3_neuronlayer.py

"""

import numpy as np
from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer

def sentence_to_ascii_list(sentence):
    ascii_list = [ord(char) for char in sentence]
    return ascii_list

def ascii_list_to_sentence(ascii_list):
    sentence = ''.join(chr(value) for value in ascii_list if value != 0)
    return sentence

def pad_list(ascii_list, max_length):
    return ascii_list + [0] * (max_length - len(ascii_list))

sentence = "Hi. How can I help you today?"+" "*47000+"\n"
ascii_list = sentence_to_ascii_list(sentence)

max_length = 50256
padded_list = pad_list(ascii_list, max_length)

input_size = max_length
target_size = max_length

ds = SupervisedDataSet(input_size, target_size)

ds.addSample(padded_list, padded_list)

net = buildNetwork(input_size, 640, 640, target_size, bias=True)

trainer = BackpropTrainer(net, ds, learningrate=0.0007515879227143)

list(map(lambda a: print(f"(Run :{a+1}) -- (Total error: {trainer.train()})"), range(60))).pop()

output = net.activate(padded_list)
processed_list = [int(abs(round(value))) for value in output]

unpadded_processed_list = [value for value in processed_list if value != 0]

processed_sentence = ascii_list_to_sentence(unpadded_processed_list)
print(processed_sentence)

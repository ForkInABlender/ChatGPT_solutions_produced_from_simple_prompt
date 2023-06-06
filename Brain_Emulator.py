# Dylan Kenneth Eliot & GPT-4


"""
After fairful consideration, after carefule consideration, my initial brain emulator setup model wasn't good
 enough in my eyes. It just needed that special loving touch of math. So I made it simpler. Now it is a 
  matter of time, compute, and data size. Well, and some eeg data about our brains using neuropy.

"""


import torch
import numpy as np
from transformers import GPT2LMHeadModel, GPT2Tokenizer
from pybrain3.structure import FeedForwardNetwork
from pybrain3.tools.shortcuts import buildNetwork

def get_gpt2_model_and_tokenizer(model_name):
    model = GPT2LMHeadModel.from_pretrained(model_name)
    tokenizer = GPT2Tokenizer.from_pretrained(model_name)
    return model, tokenizer

def create_pybrain3_network(input_size, hidden_size, output_size):
    network = buildNetwork(input_size, hidden_size, output_size, bias=True)
    return network

def gpt2_to_pybrain3(gpt2_output, pybrain3_network):
    gpt2_output_np = gpt2_output.detach().numpy()
    gpt2_output_flat = gpt2_output_np.flatten()
    pybrain3_output = pybrain3_network.activate(gpt2_output_flat)
    return pybrain3_output

def pybrain3_to_gpt2(pybrain3_output):
    pybrain3_output_torch = torch.from_numpy(pybrain3_output).float().unsqueeze(0)
    return pybrain3_output_torch

def generate_text(gpt2_model, gpt2_tokenizer, input_text, max_length=100, num_return_sequences=1):
    input_ids = gpt2_tokenizer.encode(input_text, return_tensors="pt")
    with torch.no_grad():
        output = gpt2_model.generate(input_ids, max_length=max_length, num_return_sequences=num_return_sequences)
    output_text = gpt2_tokenizer.decode(output[0], skip_special_tokens=True)
    return output_text

# Load the GPT-2 models and create the PyBrain3 network
model1, tokenizer1 = get_gpt2_model_and_tokenizer("gpt2-medium")
model2, tokenizer2 = get_gpt2_model_and_tokenizer("gpt2-large")
amygdala = create_pybrain3_network(768, 13385, 30)

# Generate output from the first GPT-2 model
input_text = "Hello, how are you?"
input_ids = tokenizer1.encode(input_text, return_tensors="pt")
with torch.no_grad():
    output1, hidden_states = model1(input_ids, output_hidden_states=True)

# Pass the GPT-2 output to PyBrain3 and back to GPT-2 format
pybrain3_output = gpt2_to_pybrain3(hidden_states[-1], amygdala)
gpt2_input = pybrain3_to_gpt2(pybrain3_output)

# Generate output from the second GPT-2 model
final_output_text = generate_text(model2, tokenizer2, gpt2_input)

print(final_output_text)

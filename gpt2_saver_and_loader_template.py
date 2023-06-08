# Dylan Kenneth Eliot & GPT-3
#
#

import torch
from transformers import GPT2Tokenizer, GPT2Model

# Initialize GPT-2 model
model_name = "gpt2"
model = GPT2Model.from_pretrained(model_name)
tokenizer = GPT2Tokenizer.from_pretrained(model_name)

# Save the GPT-2 model
output_path = "path/to/save/model.pt"
torch.save(model.state_dict(), output_path)

# Load the GPT-2 model
loaded_model = GPT2Model.from_pretrained(model_name)
loaded_model.load_state_dict(torch.load(output_path))
loaded_model.eval()

# Generate text using the loaded model
input_text = "Hello, how are you?"
inputs = tokenizer.encode(input_text, return_tensors="pt")
with torch.no_grad():
    outputs = loaded_model.generate(inputs)

# Decode and print the generated text
generated_text = tokenizer.decode(outputs[0])
print(generated_text)

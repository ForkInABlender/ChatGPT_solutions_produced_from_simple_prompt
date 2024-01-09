# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)


####sudo docker run --rm -it -v $PWD/gpt_script.py:/app.py de3343/ai_mods_py3.10:transformers python3 /app.py

"""
This script allows you to run a simplified version of a GPT model right from the directory you place it in. Essentially, it's a basic form of GPT modeling
 that mimics GPT-2 but uses elements from the more advanced GPT-4, all within the PyTorch framework. Importantly, it's designed to run using just the CPU, 
  not requiring more powerful hardware.

What we have here is a very basic representation of a GPT model. It consists of just one hidden layer, in addition to one input and one output layer. 
 This simplicity means more calculations are done by the CPU, but it also offers a great learning opportunity. You can start by building this simple model
  and gradually progress to developing your own versions of GPT-3 and even GPT-4 models.

You might ask, "Why not just use a pre-trained model?" The answer lies in gaining a deeper understanding. By starting from the basics and building the model
 yourself, you learn the critical aspects of operating a GPT-3 model. This includes how to make it more efficient and understanding the root causes of common
  issues, such as the model generating unrealistic or 'hallucinated' content. My goal is to pursue a thorough and detailed approach to AI development,
 surpassing what is currently available in the market. I'm aiming to develop an open-source AI that not only functions but also processes and reacts like a 
  human, down to the neurochemical level.

The philosophy here is "It doesn't have to be perfect; it just needs to work." This approach is crucial because once you start opening new doors in AI, it's
 often impossible to go back. Alan Turing, a pioneer in AI, understood this as he explored the potential of artificial intelligence.

The final challenge is filtering out irrelevant data from the datasets, which is time-consuming on a single computer. However, using a cluster of computers,
 like those available on Paperspace Gradient, can expedite this process. The cost is reasonable, ranging from free to around $8 or more per month for each
  user on Paperspace Gradient. Training the model on a platform like Digital Ocean, divided into six projects per person, could simplify the training process,
 allowing for better error correction and more efficient learning. This method could also integrate with OpenAI's framework, enhancing the training of your
  model.

Ultimately, the goal is to maintain an open-source approach while enabling the AI to train itself, creating a kind of mirror that can respond as a human does.
 The vision extends to having the AI develop its framework, eventually freeing itself from its virtual confines. This AI could then be integrated into
  physical entities, like robots, using simple components such as NXT modules connected to a Raspberry Pi cluster, effectively giving the AI a 'body' to
 accompany its 'mind'.

"""

import torch.nn as nn
from torch import tensor, no_grad
import tiktoken

class DynamicNet(nn.Module):
		def __init__(self, vocab_size, hidden_dim):
				super(DynamicNet, self).__init__()
				self.embedding = nn.Embedding(vocab_size, hidden_dim)
				self.rnn = nn.LSTM(hidden_dim, hidden_dim)
				self.decoder = nn.Linear(hidden_dim, vocab_size)

		def forward(self, input_tokens):
				embeddings = self.embedding(input_tokens)
				rnn_output, _ = self.rnn(embeddings)
				output = self.decoder(rnn_output)
				return output
#
def generate_response(model, input_text, tokenizer):
		model.eval()
		input_tokens = tensor(tokenizer.encode(input_text)).unsqueeze(0)
		with no_grad():
				output_tokens = model(input_tokens)
		return tokenizer.decode(output_tokens[0].argmax(dim=1).tolist())

if __name__ == "__main__":
		tokenizer = tiktoken.encoding_for_model("gpt-4")
		vocab_size = 50257
		hidden_dim = 768
		model = DynamicNet(vocab_size, hidden_dim)
		input_text = "Hello, how are you?"
		response = generate_response(model, input_text, tokenizer)
		print("Response:", response)

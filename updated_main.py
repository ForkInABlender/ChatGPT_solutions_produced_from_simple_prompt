# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""
This is similar to the gpt_script.py and functions the exact same way but with far more training needed.

This is tiktoken and pybrain3.

Tiktoken does the encoding, numpy parts do the work for BPE of the response, then it gets decoded. With no change of input, the output
 can be expected to be the same. As it also needs a feed back loop, it probably is a good idea to use recurrent or concurrent networks,
  or both. All depends on use-case.

Never-the-less this is another examle of how gpt-2 functions at its core outside its normal housing.

Because it is using lightweight mechanical modeling for static computational graphs for something dynamic, it must also pay attention to
 to its own formulation of response.

From here, it is a matter of finding a dataset that doesn't use fake or stolen data/content.
"""

from pybrain3.structure import LinearLayer, TanhLayer, LSTMLayer, SoftmaxLayer, FeedForwardNetwork, FullConnection
import tiktoken

class DynamicNet(FeedForwardNetwork):
    def __init__(self, vocab_size, hidden_dim):
        super(DynamicNet, self).__init__()
        # Layers
        self.inLayer = LinearLayer(vocab_size)
        self.tanhLayer = TanhLayer(hidden_dim)
        self.lstmLayer = LSTMLayer(hidden_dim)
        self.outLayer = LinearLayer(vocab_size)
        # Adding modules and connections
        self.addInputModule(self.inLayer)
        self.addModule(self.tanhLayer)
        self.addModule(self.lstmLayer)
        self.addOutputModule(self.outLayer)
        self.in_to_tanh = FullConnection(self.inLayer, self.tanhLayer)
        self.tanh_to_lstm = FullConnection(self.tanhLayer, self.lstmLayer)
        self.lstm_to_out = FullConnection(self.lstmLayer, self.outLayer)
        self.addConnection(self.in_to_tanh)
        self.addConnection(self.tanh_to_lstm)
        self.addConnection(self.lstm_to_out)
        self.sortModules()
#
def generate_response(model, input_text, tokenizer):
    input_tokens = tokenizer.encode(input_text)
    input_vector = [0] * vocab_size
    for token in input_tokens:
        input_vector[token] = 1
    output_tokens = model.activate(input_vector)
    return tokenizer.decode([output_tokens.argmax(axis=0)])

if __name__ == "__main__":
    tokenizer = tiktoken.encoding_for_model("gpt-4")
    vocab_size = 50257
    hidden_dim = 768
    model = DynamicNet(vocab_size, hidden_dim)
    input_text = "Hello, how are you?"
    response = generate_response(model, input_text, tokenizer)
    print("Response:", response)

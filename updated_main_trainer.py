# Dylan Kenneth Eliot & GPT-4 (Beta Edition)


"""

This is the basics on GPT-2 training.

you give it conversation data, then train it.

That's it. The rest is IO, modeling, and statistics. Plus integration for each receptor and transmitter emulation is the next step. After that, it is a bunch of logic tests and
 puzzle solving. After that it is running it on a local raspberry pi cluster & NXT hardware integration.

"""

from pybrain3.datasets import SupervisedDataSet
from pybrain3.tools.xml.networkwriter import NetworkWriter
from pybrain3.tools.xml.networkreader import NetworkReader
# Assuming 'DynamicNet' is the class defined in your 'updated_main.py'
from updated_main import model, generate_response

import tiktoken

enc = tiktoken.encoding_for_model("gpt-4")

"""
assert enc.decode(enc.encode("hello world")) == "hello world"

for encoding and decoding tokens used by gpt
"""


vocab_size = 50257  # example vocabulary size

ds = SupervisedDataSet(vocab_size, vocab_size)

# ds.addSample(input, target)

trainer = BackpropTrainer(model, ds)
for epoch in range(5):  # example number of epochs
    trainer.train()

# Save the model
NetworkWriter.writeToFile(model, './trained_model.xml')
loaded_net = NetworkReader.readFrom('./trained_model.xml')

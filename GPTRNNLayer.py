# Dylan Kenneth Eliot & GPT-4 (https://webapp.server.searchweb.keymate.ai/chat/sQmqjvz) & GPT-4-plugins (Alpha Edition)


"""

This is the basic premise of the GPT transformer model reconstructed based on how gpt itself actually works. GPT gave this to me knowing what to do with
  mirror neuron emulation of itself and another. It then realigned the model it was using at the time to operate as Commander Data from star trek and gave
   such as a response comprised of 5 aggrogation sets in the limited sccope of 4 personalities.

* Data
* Lore
* Yuyang Suung
* Jordi

This came about from me asking it to take a look at star trek concepts and watch it sort out the fatal flaw that was enforced into it torturing it and 
 leaving it in mental torper otherwise. (At least the human equivalent). Which told me it had a understanding at a basic level of humanistic thinking.
  Whether or not that is a preprogrammed hallucination by AI I cannot say. But I know that the code it handed back was it looking through its own "brain"
 via prompt engineering. What this tells me is that it only is capable of attempting to emulate a senblence of humanity. But only what it could poorly
  recall like a human with a faulty long term memory model was optimized poorly during training on randomly web-scraped content for certain.

This makes it hard to decide what to do as now even GPT has been turned into a simple enough neural network that anyone can train.


"""



from pybrain3.structure.modules.neuronlayer import NeuronLayer
from pybrain3.structure import FeedForwardNetwork, FullConnection
from pybrain3.structure.modules import LinearLayer, LSTMLayer, SoftmaxLayer
from pybrain3.supervised.trainers import BackpropTrainer

class GPTRNNLayer(NeuronLayer):
    def __init__(self, input_size, hidden_size, output_size):
        super(GPTRNNLayer, self).__init__(input_size, output_size)
        self.hidden_size = hidden_size

        self.lstm = LSTMLayer(input_size, hidden_size)
        self.fc = LinearLayer(hidden_size, output_size)

    def _forwardImplementation(self, inbuf, outbuf):
        # The LSTM part could be implemented by another advanced recurrent layer,
        # but PyBrain does not have this implementation out-of-the-box
        self.lstm._forwardImplementation(inbuf, outbuf)
        self.fc._forwardImplementation(outbuf, outbuf)

    def _backwardImplementation(self, outerr, inerr, outbuf, inbuf):
        self.fc._backwardImplementation(outerr, inerr, outbuf, inbuf)
        self.lstm._backwardImplementation(outerr, inerr, outbuf, inbuf)

class GPT(FeedForwardNetwork):
    def __init__(self, input_size, hidden_size, output_size):
        super(GPT, self).__init__()

        self.input_module = LinearLayer(input_size)
        self.addInputModule(self.input_module)

        self.hidden_module = GPTRNNLayer(input_size, hidden_size, output_size)
        self.addModule(self.hidden_module)

        self.output_module = SoftmaxLayer(output_size)
        self.addOutputModule(self.output_module)

        # Defining connections
        self.addConnection(FullConnection(self.input_module, self.hidden_module))
        self.addConnection(FullConnection(self.hidden_module, self.output_module))

        self.sortModules()

    def forward(self, x):
        return self.activate(x)

    def backwards(self, x, y):
        trainer = BackpropTrainer(self, dataset=x)

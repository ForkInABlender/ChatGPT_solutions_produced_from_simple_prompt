# Dylan Kenneth Eliot
"""
This creates with the Dim3NeuronLayer as the hidden layers. This should form the basic template for GPT knowing that the Dim3NeuronLayer does the QKV as well as 3d, 2d, & 1d processing of data.

The Dim3NeuronLayer itself acts as a layer for allowing a 3 dimensional neuron comprised of its 2 dimensional counter parts, using them as anchor corner points.
 This gives the FFT the ability to update the connections and individual data points, allowing for 3d processing, spatial awareness, image processing, and language processing.

language datasets soon to follow as well as brain-emulator code revision.

"""

from pybrain3.structure import FeedForwardNetwork, LinearLayer, BiasUnit
from pybrain3.structure.connections import FullConnection
##############################################
from dim3_neuronlayer import Dim3NeuronLayer #
##############################################
#from pybrain3.datasets import SupervisedDataSet

net = FeedForwardNetwork()

input_size = 50257
output_size = 50257
num_hidden_layers = 30 # Max number on (768, 768, 768) spatial layout per 3d hidden layer

indim, outdim, num_heads = 768, 768, 256

in_layer = LinearLayer(input_size)
out_layer = LinearLayer(output_size)
bias = BiasUnit()

net.addInputModule(in_layer)
net.addOutputModule(out_layer)
net.addModule(bias)

previous_layer = in_layer

for i in range(num_hidden_layers):
    hidden_layer = Dim3NeuronLayer(indim=indim, outdim=outdim, num_heads=num_heads, name=f"HiddenLayer_{i+1}")
    net.addModule(hidden_layer)

    net.addConnection(FullConnection(previous_layer, hidden_layer))
    net.addConnection(FullConnection(bias, hidden_layer))

    previous_layer = hidden_layer

net.addConnection(FullConnection(previous_layer, out_layer))
net.addConnection(FullConnection(bias, out_layer))

net.sortModules()

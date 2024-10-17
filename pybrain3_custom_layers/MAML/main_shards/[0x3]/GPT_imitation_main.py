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


if __name__ == "__main__":
    # Create a feed-forward network
    net = FeedForwardNetwork()
    dim_s=768
    layers=96

    # Input layer
    in_layer = LinearLayer(50257)
    net.addInputModule(in_layer)

    hidden_layers = []
    for i in range(layers):
        hidden_layer = Dim3NeuronLayer(dim_s, dim_s, embed_dropout=0, attn_dropout=0, activation_function='gelu', lr=0.0002, weight_decay=0.1, gradient_clipping=1)
        net.addModule(hidden_layer)
        hidden_layers.append(hidden_layer)

    # Output layer
    out_layer = LinearLayer(50257)
    net.addOutputModule(out_layer)

    # Full connections
    in_to_hidden = FullConnection(in_layer, hidden_layers[0])
    net.addConnection(in_to_hidden)

    for i in range(layers-1):
        net.addConnection(FullConnection(hidden_layers[i], hidden_layers[i + 1]))

    hidden_to_out = FullConnection(hidden_layers[-1], out_layer)
    net.addConnection(hidden_to_out)

    # Finalize the network
    net.sortModules()
    ############

    # Example input
    input_data = [0.001] * 50257
    output_data = net.activate(input_data)
    print("Output:", output_data)

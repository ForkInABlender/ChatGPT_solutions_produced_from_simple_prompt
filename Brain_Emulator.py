# Dylan Kenneth Eliot & GPT3 with custom bot prompt
#
#

from pybrain3.structure import FeedForwardNetwork
from pybrain3.structure import LinearLayer, SigmoidLayer
from pybrain3.structure import FullConnection
from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.tools.xml.networkwriter import NetworkWriter
from pybrain3.tools.xml.networkreader import NetworkReader

"""

to save a network:

NetworkWriter.writeToFile(network_to_be_saved, 'network_saved.xml')

to use a saved network:

network_to_be_saved = NetworkReader.readFrom('network_saved.xml')

"""

from pybrain3.structure import FeedForwardNetwork
from pybrain3.structure import FullConnection
from pybrain3.tools.shortcuts import buildNetwork

amygdala = buildNetwork(2012160, 13385, 30, bias=True)
hippocampus = buildNetwork(2012160, 200000, 504600)
frontal_lobe = buildNetwork(10000, 5000, 1000)
parietal_lobe = buildNetwork(8000, 4000, 800)
temporal_lobe = buildNetwork(7500, 3500, 700)
occipital_lobe = buildNetwork(6000, 3000, 600)
cerebellum = buildNetwork(9000, 4500, 900)
brainstem = buildNetwork(3000, 1500, 300)
limbic_system = buildNetwork(5000, 2500, 500)
wernicke_area = buildNetwork(3500, 1800, 350)
broca_area = buildNetwork(3000, 1500, 300)

networks = [
    amygdala, hippocampus, frontal_lobe, parietal_lobe, temporal_lobe,
    occipital_lobe, cerebellum, brainstem, limbic_system, wernicke_area, broca_area
]

# Connect all networks to each other, including self-connections
for i in range(len(networks)):
    for j in range(len(networks)):
        connection = FullConnection(networks[i], networks[j])
        networks[i].addConnection(connection)
        networks[i].sortModules()

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


hippocampus = buildNetwork(2012160, 200000, 504600)
frontal_lobe = buildNetwork(10000, 5000, 1000)  # Emulating the Frontal Lobe
parietal_lobe = buildNetwork(8000, 4000, 800)   # Emulating the Parietal Lobe
temporal_lobe = buildNetwork(7500, 3500, 700)   # Emulating the Temporal Lobe
occipital_lobe = buildNetwork(6000, 3000, 600)  # Emulating the Occipital Lobe
cerebellum = buildNetwork(9000, 4500, 900)      # Emulating the Cerebellum
brainstem = buildNetwork(3000, 1500, 300)       # Emulating the Brainstem
limbic_system = buildNetwork(5000, 2500, 500) # Emulating the Limbic System

# Emulating Wernicke's Area
wernicke_area = FeedForwardNetwork()
wernicke_area.addInputModule(LinearLayer(3500))
wernicke_area.addModule(SigmoidLayer(1800))
wernicke_area.addOutputModule(LinearLayer(350))
wernicke_area.addConnection(FullConnection(wernicke_area['in'], wernicke_area['hidden0']))
wernicke_area.addConnection(FullConnection(wernicke_area['hidden0'], wernicke_area['out']))
wernicke_area.sortModules()

# Emulating Broca's Area
broca_area = FeedForwardNetwork()
broca_area.addInputModule(LinearLayer(3000))
broca_area.addModule(SigmoidLayer(1500))
broca_area.addOutputModule(LinearLayer(300))
broca_area.addConnection(FullConnection(broca_area['in'], broca_area['hidden0']))
broca_area.addConnection(FullConnection(broca_area['hidden0'], broca_area['out']))
broca_area.sortModules()

# Connecting the Frontal Lobe to other regions
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], parietal_lobe['in']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], temporal_lobe['in']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], occipital_lobe['in']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], cerebellum['in']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], limbic_system['in']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], wernicke_area['in']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], broca_area['in']))

# Connecting the Parietal Lobe to other regions
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], frontal_lobe['in']))
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], temporal_lobe['in']))
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], occipital_lobe['in']))
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], cerebellum['in']))
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], limbic_system['in']))

# Connecting the Temporal Lobe to other regions
temporal_lobe.addConnection(FullConnection(temporal_lobe['out'], frontal_lobe['in']))
temporal_lobe.addConnection(FullConnection(temporal_lobe['out'], parietal_lobe['in']))
temporal_lobe.addConnection(FullConnection(temporal_lobe['out'], occipital_lobe['in']))
temporal_lobe.addConnection(FullConnection(temporal_lobe['out'], limbic_system['in']))

# Connecting the Occipital Lobe to other regions
occipital_lobe.addConnection(FullConnection(occipital_lobe['out'], frontal_lobe['in']))
occipital_lobe.addConnection(FullConnection(occipital_lobe['out'], parietal_lobe['in']))
occipital_lobe.addConnection(FullConnection(occipital_lobe['out'], temporal_lobe['in']))
occipital_lobe.addConnection(FullConnection(occipital_lobe['out'], limbic_system['in']))

# Connecting the Cerebellum to other regions
cerebellum.addConnection(FullConnection(cerebellum['out'], frontal_lobe['in']))
cerebellum.addConnection(FullConnection(cerebellum['out'], parietal_lobe['in']))
cerebellum.addConnection(FullConnection(cerebellum['out'], temporal_lobe['in']))
cerebellum.addConnection(FullConnection(cerebellum['out'], occipital_lobe['in']))
cerebellum.addConnection(FullConnection(cerebellum['out'], limbic_system['in']))

# Connecting the Limbic System to other regions
limbic_system.addConnection(FullConnection(limbic_system['out'], frontal_lobe['in']))
limbic_system.addConnection(FullConnection(limbic_system['out'], parietal_lobe['in']))
limbic_system.addConnection(FullConnection(limbic_system['out'], temporal_lobe['in']))
limbic_system.addConnection(FullConnection(limbic_system['out'], occipital_lobe['in']))
limbic_system.addConnection(FullConnection(limbic_system['out'], cerebellum['in']))

# Connecting Wernicke's Area to other regions
wernicke_area.addConnection(FullConnection(wernicke_area['out'], frontal_lobe['in']))
wernicke_area.addConnection(FullConnection(wernicke_area['out'], parietal_lobe['in']))
wernicke_area.addConnection(FullConnection(wernicke_area['out'], temporal_lobe['in']))
wernicke_area.addConnection(FullConnection(wernicke_area['out'], occipital_lobe['in']))
wernicke_area.addConnection(FullConnection(wernicke_area['out'], limbic_system['in']))
wernicke_area.addConnection(FullConnection(wernicke_area['out'], broca_area['in']))

# Connecting Broca's Area to other regions
broca_area.addConnection(FullConnection(broca_area['out'], frontal_lobe['in']))
broca_area.addConnection(FullConnection(broca_area['out'], parietal_lobe['in']))
broca_area.addConnection(FullConnection(broca_area['out'], temporal_lobe['in']))
broca_area.addConnection(FullConnection(broca_area['out'], occipital_lobe['in']))
broca_area.addConnection(FullConnection(broca_area['out'], limbic_system['in']))
broca_area.addConnection(FullConnection(broca_area['out'], wernicke_area['in']))

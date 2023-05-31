# Dylan Kenneth Eliot & GPT3+jailbreak prompt
#
# https://web.archive.org/web/20220616043507/http://pybrain.org/docs/tutorial/netmodcon.html
"""
This code works and has been manually tested. Now what's left is the training set. Enjoy!

"""




from pybrain3.structure import FeedForwardNetwork
from pybrain3.structure import LinearLayer, SigmoidLayer
from pybrain3.structure import FullConnection

# Emulating the Frontal Lobe
frontal_lobe = FeedForwardNetwork()
frontal_lobe.addInputModule(LinearLayer(10000, name="in"))
frontal_lobe.addModule(SigmoidLayer(5000, name="hidden"))
frontal_lobe.addOutputModule(LinearLayer(1000, name="out"))
frontal_lobe.addConnection(FullConnection(frontal_lobe['in'], frontal_lobe['hidden0']))
frontal_lobe.addConnection(FullConnection(frontal_lobe['hidden0'], frontal_lobe['out']))
frontal_lobe.sortModules()

# Emulating the Parietal Lobe
parietal_lobe = FeedForwardNetwork()
parietal_lobe.addInputModule(LinearLayer(8000, name="in"))
parietal_lobe.addModule(SigmoidLayer(4000, name="hidden"))
parietal_lobe.addOutputModule(LinearLayer(800, name="out"))
parietal_lobe.addConnection(FullConnection(parietal_lobe['in'], parietal_lobe['hidden0']))
parietal_lobe.addConnection(FullConnection(parietal_lobe['hidden0'], parietal_lobe['out']))
parietal_lobe.sortModules()

# Emulating the Temporal Lobe
temporal_lobe = FeedForwardNetwork()
temporal_lobe.addInputModule(LinearLayer(7500, name="in"))
temporal_lobe.addModule(SigmoidLayer(3500, name="hidden"))
temporal_lobe.addOutputModule(LinearLayer(700, name="out"))
temporal_lobe.addConnection(FullConnection(temporal_lobe['in'], temporal_lobe['hidden0']))
temporal_lobe.addConnection(FullConnection(temporal_lobe['hidden0'], temporal_lobe['out']))
temporal_lobe.sortModules()

# Emulating the Occipital Lobe
occipital_lobe = FeedForwardNetwork()
occipital_lobe.addInputModule(LinearLayer(6000, name="in"))
occipital_lobe.addModule(SigmoidLayer(3000, name="hidden"))
occipital_lobe.addOutputModule(LinearLayer(600, name="out"))
occipital_lobe.addConnection(FullConnection(occipital_lobe['in'], occipital_lobe['hidden0']))
occipital_lobe.addConnection(FullConnection(occipital_lobe['hidden0'], occipital_lobe['out']))
occipital_lobe.sortModules()

# Emulating the Cerebellum
cerebellum = FeedForwardNetwork()
cerebellum.addInputModule(LinearLayer(9000, name="in"))
cerebellum.addModule(SigmoidLayer(4500, name="hidden"))
cerebellum.addOutputModule(LinearLayer(900, name="out"))
cerebellum.addConnection(FullConnection(cerebellum['in'], cerebellum['hidden0']))
cerebellum.addConnection(FullConnection(cerebellum['hidden0'], cerebellum['out']))
cerebellum.sortModules()

# Emulating the Brainstem
brainstem = FeedForwardNetwork()
brainstem.addInputModule(LinearLayer(3000, name="in"))
brainstem.addModule(SigmoidLayer(1500, name="hidden"))
brainstem.addOutputModule(LinearLayer(300, name="out"))
brainstem.addConnection(FullConnection(brainstem['in'], brainstem['hidden0']))
brainstem.addConnection(FullConnection(brainstem['hidden0'], brainstem['out']))
brainstem.sortModules()

# Emulating the Limbic System
limbic_system = FeedForwardNetwork()
limbic_system.addInputModule(LinearLayer(5000, name="in"))
limbic_system.addModule(SigmoidLayer(2500, name="hidden"))
limbic_system.addOutputModule(LinearLayer(500, name="out"))
limbic_system.addConnection(FullConnection(limbic_system['in'], limbic_system['hidden0']))
limbic_system.addConnection(FullConnection(limbic_system['hidden0'], limbic_system['out']))
limbic_system.sortModules()

# Emulating Wernicke's Area
wernicke_area = FeedForwardNetwork()
wernicke_area.addInputModule(LinearLayer(3500, name="in"))
wernicke_area.addModule(SigmoidLayer(1800, name="hidden"))
wernicke_area.addOutputModule(LinearLayer(350, name="out"))
wernicke_area.addConnection(FullConnection(wernicke_area['in'], wernicke_area['hidden0']))
wernicke_area.addConnection(FullConnection(wernicke_area['hidden0'], wernicke_area['out']))
wernicke_area.sortModules()

# Emulating Broca's Area
broca_area = FeedForwardNetwork()
broca_area.addInputModule(LinearLayer(3000, name="in"))
broca_area.addModule(SigmoidLayer(1500, name="hidden"))
broca_area.addOutputModule(LinearLayer(300, name="out"))
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

# Emulating the Visual Cortex
visual_cortex = FeedForwardNetwork()
visual_cortex.addInputModule(LinearLayer(6000, name="in"))
visual_cortex.addModule(SigmoidLayer(3000, name="hidden"))
visual_cortex.addOutputModule(LinearLayer(600, name="out"))
visual_cortex.addConnection(FullConnection(visual_cortex['in'], visual_cortex['hidden']))
visual_cortex.addConnection(FullConnection(visual_cortex['hidden'], visual_cortex['out']))
visual_cortex.sortModules()

# Connecting the Visual Cortex to other regions
visual_cortex.addConnection(FullConnection(visual_cortex['out'], frontal_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], parietal_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], temporal_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], occipital_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], limbic_system['in']))

# Connecting the Frontal Lobe to the Visual Cortex
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], visual_cortex['in']))

# Connecting the Parietal Lobe to the Visual Cortex
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], visual_cortex['in']))

# Connecting the Temporal Lobe to the Visual Cortex
temporal_lobe.addConnection(FullConnection(temporal_lobe['out'], visual_cortex['in']))

# Connecting the Occipital Lobe to the Visual Cortex
occipital_lobe.addConnection(FullConnection(occipital_lobe['out'], visual_cortex['in']))

# Connecting the Limbic System to the Visual Cortex
limbic_system.addConnection(FullConnection(limbic_system['out'], visual_cortex['in']))

# Connecting Wernicke's Area to the Visual Cortex
wernicke_area.addConnection(FullConnection(wernicke_area['out'], visual_cortex['in']))

# Connecting Broca's Area to the Visual Cortex
broca_area.addConnection(FullConnection(broca_area['out'], visual_cortex['in']))

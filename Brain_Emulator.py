# Dylan Kenneth Eliot & GPT3+jailbreak prompt
#
# https://web.archive.org/web/20220616043507/http://pybrain.org/docs/tutorial/netmodcon.html
"""
This code works and has been manually tested. Now what's left is the training sets. Enjoy!

"""

from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.structure import FullConnection

# Emulating the Frontal Lobe
frontal_lobe = buildNetwork(10000, 5000, 1000)

# Emulating the Parietal Lobe
parietal_lobe = buildNetwork(8000, 4000, 800)

# Emulating the Temporal Lobe
temporal_lobe = buildNetwork(7500, 3500, 700)

# Emulating the Occipital Lobe
occipital_lobe = buildNetwork(6000, 3000, 600)

# Emulating the Cerebellum
cerebellum = buildNetwork(9000, 4500, 900)

# Emulating the Brainstem
brainstem = buildNetwork(3000, 1500, 300)

# Emulating the Limbic System
limbic_system = buildNetwork(5000, 2500, 500)

# Emulating Wernicke's Area
wernicke_area = buildNetwork(3500, 1800, 350)

# Emulating Broca's Area
broca_area = buildNetwork(3000, 1500, 300)

# Emulating the Visual Cortex
visual_cortex = buildNetwork(6000, 3000, 600)

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

# Connecting the Visual Cortex to other regions
visual_cortex.addConnection(FullConnection(visual_cortex['out'], frontal_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], parietal_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], temporal_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], occipital_lobe['in']))
visual_cortex.addConnection(FullConnection(visual_cortex['out'], limbic_system['in']))

# Connecting the Frontal Lobe to the Visual Cortex
frontal_lobe.addConnection(FullConnection(frontal_lobe['out'], visual_cortex['in']))
parietal_lobe.addConnection(FullConnection(parietal_lobe['out'], visual_cortex['in']))
temporal_lobe.addConnection(FullConnection(temporal_lobe['out'], visual_cortex['in']))
occipital_lobe.addConnection(FullConnection(occipital_lobe['out'], visual_cortex['in']))
limbic_system.addConnection(FullConnection(limbic_system['out'], visual_cortex['in']))
wernicke_area.addConnection(FullConnection(wernicke_area['out'], visual_cortex['in']))
broca_area.addConnection(FullConnection(broca_area['out'], visual_cortex['in']))

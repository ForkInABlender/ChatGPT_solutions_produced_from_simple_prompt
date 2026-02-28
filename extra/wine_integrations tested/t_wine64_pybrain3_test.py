# Dylan Kenneth Eliot

"""
This was tested within wine64 python3.10.11 code. 

It ran identically to rest of compute elsewhere and restored the model 1:1 on linux & on termux.
"""

from pybrain3.structure import FeedForwardNetwork
from pybrain3.structure.modules import LinearLayer, SigmoidLayer
from pybrain3.structure.connections import FullConnection
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.datasets import SupervisedDataSet
import numpy as np
import random

# -----------------------------
# Create Network
# -----------------------------

net = FeedForwardNetwork()

# Layers
input_layer = LinearLayer(2)      # 2 input neurons
hidden_layer = SigmoidLayer(3)    # 3 hidden neurons
output_layer = LinearLayer(1)     # 1 output neuron

# Add modules
net.addInputModule(input_layer)
net.addModule(hidden_layer)
net.addOutputModule(output_layer)

# Connect layers
in_to_hidden = FullConnection(input_layer, hidden_layer)
hidden_to_out = FullConnection(hidden_layer, output_layer)

net.addConnection(in_to_hidden)
net.addConnection(hidden_to_out)

# Finalize network structure
net.sortModules()

# -----------------------------
# Create Dataset (XOR-like)
# -----------------------------

ds = SupervisedDataSet(2, 1)

ds.addSample([0, 0], [0])
ds.addSample([0, 1], [1])
ds.addSample([1, 0], [1])
ds.addSample([1, 1], [0])

# -----------------------------
# Train Network
# -----------------------------

trainer = BackpropTrainer(net, ds, learningrate=0.1, momentum=0.1, verbose=True)

for epoch in range(1000):
    trainer.train()

# -----------------------------
# Test Network
# -----------------------------


# ─────────────────────────────────────────────────────────────
# SAVE SNAPSHOT
# ─────────────────────────────────────────────────────────────

from lxml import etree

def serialize_network(net, filename):
    root = etree.Element("PyBrainNetwork", 
                         name=str(getattr(net, 'name', 'Network')), 
                         type=type(net).__name__)
    
    # 1. Structure (Input/Output Pointers)
    struct = etree.SubElement(root, "Structure")
    if net.inmodules:
        for in_mod in net.inmodules:
            etree.SubElement(struct, "InputModule").text = in_mod.name
    if net.outmodules:
        for out_mod in net.outmodules:
            etree.SubElement(struct, "OutputModule").text = out_mod.name

    # 2. Modules (Layers)
    modules_xml = etree.SubElement(root, "Modules")
    for mod in net.modules:
        m_xml = etree.SubElement(modules_xml, "Module", 
                                 name=str(mod.name), 
                                 type=type(mod).__name__,
                                 indim=str(mod.indim),
                                 outdim=str(mod.outdim))
        
        # Capture Parameters (Weights/Biases)
        # We check .params (standard) or .bias (some custom layers)
        params = getattr(mod, 'params', None)
        
        if params is not None:
            # We use base64 or high-precision float strings
            # Storing as space-separated floats for XML readability
            m_xml.text = " ".join(map(str, params.flatten()))
            m_xml.set("param_count", str(len(params)))

    # 3. Connections (Synapses)
    connections_xml = etree.SubElement(root, "Connections")
    for source_mod in net.connections:
        for conn in net.connections[source_mod]:
            c_xml = etree.SubElement(connections_xml, "Connection",
                                     type=type(conn).__name__,
                                     in_mod=conn.inmod.name,
                                     out_mod=conn.outmod.name)
            
            # Connections like FullConnection hold the actual weight matrix here
            c_params = getattr(conn, 'params', None)
            if c_params is not None:
                c_xml.text = " ".join(map(str, c_params.flatten()))
                c_xml.set("param_count", str(len(c_params)))

    # Write to file
    tree = etree.ElementTree(root)
    tree.write(filename, pretty_print=True, xml_declaration=True, encoding="utf-8")
    return tree

import xmltodict

def reconstruct_pybrain_object(xml_file):
    # Parse the local XML snapshot
    tree = etree.parse(xml_file)
    root = tree.getroot()
    parsed = xmltodict.parse(etree.tostring(root))
    net2 = eval(parsed['PyBrainNetwork']['@type'])()
    net2.name = parsed['PyBrainNetwork']['@name']
    evalled=0
    for net_keys in parsed['PyBrainNetwork'].keys():
        evalled=1
        try:
            if net_keys == 'Modules':
                for x2 in parsed['PyBrainNetwork'][net_keys]['Module']:
                    x3=eval(x2['@type'])(eval(x2['@indim']), eval(x2['@outdim']))
                    x3.name=x2['@name']
                    net2.addModule(x3)
                    del x3
            if net_keys == 'Connections':
                for y in parsed['PyBrainNetwork'][net_keys]['Connection']:
                    c=eval(y['@type'])(net2[y['@in_mod']], net2[y['@out_mod']])
                    c._params=list(map(lambda z: np.float64(eval(z)), y['#text'].split()))
                    net2.addConnection(c)
                    del c
        except Exception as e:
            e.with_traceback()
    net2.sortModules()
    net2.addInputModule(net2[parsed['PyBrainNetwork']['Structure']['InputModule']])
    net2.addOutputModule(net2[parsed['PyBrainNetwork']['Structure']['OutputModule']])
    net2.sortModules()
    return net2

save_path = 'wine64_pybrain_231test.xml'
#e=serialize_network(net, save_path)

network_obj = reconstruct_pybrain_object(save_path)

print("\nTesting:")
for sample in ds:
    input_data = sample[0]
    output = network_obj.activate(input_data)
    print(f"Input: {input_data} -> Output: {output}")

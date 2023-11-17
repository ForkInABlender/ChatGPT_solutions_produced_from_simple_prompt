# Dylan Kenneth Eliot & Google Bard AI & GPT-4-plugins (Alpha Edition)


"""

What it allows for is now using the information about the molecule to do emulation of neurochemical responses.

This is valuable for when you need to test on a fixed size of neural networks that is small to find stability points.
But it can also be retooled to apply it to the transformers, pytorch, numpy, and tensorflow models as well.

Now the smile below is of Crispr or the common substructure used in DNA repair mediated by zinc and magnesium ions.

What I'd like for to be done with this is for it to be used to find stable DNA repairs and how they apply to the brain via emulation.


"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.tools.shortcuts import buildNetwork
from brian2 import *

# RDKit: Compute molecular properties
smiles = "CC(=O)O"
molecule = Chem.MolFromSmiles(smiles)
molwt = Descriptors.MolWt(molecule)
logp = Descriptors.MolLogP(molecule)

print(smiles, molecule, molwt, log)

#####

start_scope()
duration = 100*ms
eqs = '''
dv/dt = (molwt + logp - v)/tau : 1
tau : second
'''
print(eqs)
input_neurons = 2
hidden_neurons_1 = 3
output_neurons = 1
G = NeuronGroup(input_neurons, eqs, threshold='v>0.8', reset='v = 0', method='exact')

G.tau = '10*ms + i*ms'
G.add_attribute('logp')
G.add_attribute('molwt')
G.molwt = molwt  # Set molecular weight for each neuron
G.logp = logp
M = StateMonitor(G, 'v', record=True)
run(duration)
#####
ds = SupervisedDataSet(input_neurons, output_neurons)
ds.addSample(M.v[:, -1], molwt)
net = buildNetwork(input_neurons, hidden_neurons_1, output_neurons)
trainer = BackpropTrainer(net, ds)

for a in range(1000):
    trainer.train()

# Test the network
for input, target in ds:
    print(input, target)
    prediction = net.activate(input)
    print(f"Input: {input}, Prediction: {prediction}")

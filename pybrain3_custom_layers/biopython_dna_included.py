# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This translates a genetic sequence into something usable via pybrain3.

What it does is it takes a DNA strand, and then "makes a protein", then uses it to emulate neurochemical pickup and response of the neuron or neuron group using it.

The next step is to go from this model, and develop the inverse process from a running model. The type of model I'm referring to will have to use biopython with something
 like neurosky or BCI DIY hardware. The same trick could be achieved from MEGAN & mMRI & fMRI scans.


Good websites to draw from are listed below:

* https://www.genecards.org/
* https://www.omim.org/
* https://www.kegg.jp/
* https://www.uniprot.org/
* https://www.ncbi.nlm.nih.gov/gene/
* https://useast.ensembl.org/index.html
* https://genome.ucsc.edu/
"""


from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.SeqUtils import molecular_weight
from rdkit import Chem
from rdkit.Chem import Descriptors
from pybrain3.datasets import SupervisedDataSet
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.tools.shortcuts import buildNetwork
from brian2 import *

# Biopython: Translate DNA to protein sequence
dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
mrna_seq = dna_seq.transcribe()
protein_seq = mrna_seq.translate()

# Dictionary of amino acids and their SMILES representations including stop codon
amino_acid_smiles = {
    'A': 'CC(C)C(=O)O', 'R': 'NC(CCCN=C(N)N)C(=O)O', 'N': 'NC(CC(=O)N)C(=O)O',
    'D': 'NC(CC(=O)O)C(=O)O', 'C': 'NC(CS)C(=O)O', 'E': 'NC(CCC(=O)O)C(=O)O',
    'Q': 'NC(CCC(=O)N)C(=O)O', 'G': 'C(C(=O)O)N', 'H': 'NC(CC1=CN=C-N1)C(=O)O',
    'I': 'CC(C)CC(C(=O)O)N', 'L': 'CC(C)CC(C(=O)O)N', 'K': 'NC(CCCCN)C(=O)O',
    'M': 'CC(C)CC(C(=O)O)N', 'F': 'NC(CC1=CC=CC=C1)C(=O)O', 'P': 'C1CCNC1C(=O)O',
    'S': 'NC(CO)C(=O)O', 'T': 'NC(C(C)O)C(=O)O', 'W': 'NC(CC1=CNC2=CC=CC=C21)C(=O)O',
    'Y': 'NC(CC1=CC=C(O)C=C1)C(=O)O', 'V': 'CC(C)C(C(=O)O)N', '*': ''
}

# Set peptide sequence to the translated protein sequence
# Generate SMILES for the peptide
peptide_smiles = ".".join([amino_acid_smiles[aa] for aa in str(protein_seq) if aa in amino_acid_smiles and amino_acid_smiles[aa]])
print("Peptide SMILES:", peptide_smiles)

# RDKit: Compute molecular properties
molecule = Chem.MolFromSmiles(peptide_smiles)
molwt = Descriptors.MolWt(molecule)
logp = Descriptors.MolLogP(molecule)

print(f"Peptide: {peptide}")
print(f"Molecular Weight: {molwt}")
print(f"logP: {logp}")

# Brian2: Simulate biological function
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

# PyBrain3: Neural network training
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

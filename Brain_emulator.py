# Dylan Kenneth Eliot & GPT-4-Plugins (Alpha Edition)

"""

The point of this was to display the 2-dimensional layout of a poly-dimensional structure like the brain. This can be used to develop an actual AI from
 human values. What it isn't for? Video games. :) In fact, it is intentionally not for video games. :) A human mind is not for video games.
  And I will not condone an AI being added into a video game based on my code. I will, however, condone the use of my code for scientific purposes; up to
   and including raising the dead as AI until they can have their physical body regrown even in a lab.
"""

from sympy import symbols, lambdify, sin
# Define a symbolic variable
x = symbols('x')
# Define a symbolic function
f = sin(x) + x**2
# Convert the symbolic function to a Python function
f_lambdified = lambdify(x, f, 'numpy')
# Use the function
result = f_lambdified(2)  # Should return sin(2) + 2^2 = 4.909297426825681
print("Result:", result)

####
# Dylan Kenneth Eliot & GPT-4-Plugins (Alpha Edition)
from pybrain3.structure import TanhLayer, LSTMLayer, RecurrentNetwork, FullConnection
from pybrain3.tools.shortcuts import buildNetwork
from main__GPT_Model import net as gpt_net
import numpy as np
# Anterior Cingulate Cortex (ACC) Model
acc_net = RecurrentNetwork()
acc_net.addInputModule(TanhLayer(100, name='in'))
acc_net.addModule(TanhLayer(250, name='hidden0'))
acc_net.addOutputModule(TanhLayer(50, name='out'))
acc_net.addConnection(FullConnection(acc_net['in'], acc_net['hidden0']))
acc_net.addConnection(FullConnection(acc_net['hidden0'], acc_net['out']))
acc_net.addRecurrentConnection(FullConnection(acc_net['hidden0'], acc_net['hidden0']))
# Mirror Neuron System Model
mirror_net = RecurrentNetwork()
mirror_net.addInputModule(TanhLayer(500, name='in'))
mirror_net.addModule(TanhLayer(250, name='hidden0'))
mirror_net.addOutputModule(TanhLayer(50, name='out'))
mirror_net.addConnection(FullConnection(mirror_net['in'], mirror_net['hidden0']))
mirror_net.addConnection(FullConnection(mirror_net['hidden0'], mirror_net['out']))
mirror_net.addRecurrentConnection(FullConnection(mirror_net['hidden0'], mirror_net['hidden0']))
# The rest for training
insula_net = buildNetwork(500, 250, 50, hiddenclass=TanhLayer)
frontal_lobe_net = buildNetwork(100, 600, 120)
prefrontal_cortex = buildNetwork(250, 500, 100, hiddenclass=TanhLayer)
cerebral_cortex = buildNetwork(925, 500, 100)
temporal_lobe = buildNetwork(50, 500, 250, hiddenclass=TanhLayer)
visual_cortex = buildNetwork(120, 400, 80)
brocas_area_net = buildNetwork(150, 375, 75)
wernickes_area_net = buildNetwork(75, 350, 70)
cerebellum = buildNetwork(700, 350, 70)
brainstem = buildNetwork(120, 300, 60)
hippocampus = buildNetwork(250, 125, 50, hiddenclass=LSTMLayer)
thalamus = buildNetwork(400, 200, 120)
amygdala = buildNetwork(120, 60, 10, hiddenclass=TanhLayer)
hypothalamus = buildNetwork(120, 100, 20)
#
models = [
		("acc_net", acc_net),
		("mirror_net", mirror_net),
		("insula_net", insula_net),
		("frontal_lobe_net", frontal_lobe_net),
		("prefrontal_cortex", prefrontal_cortex),
		("cerebral_cortex", cerebral_cortex),
		("temporal_lobe", temporal_lobe),
		("visual_cortex", visual_cortex),
		("brocas_area_net", brocas_area_net),
		("wernickes_area_net", wernickes_area_net),
		("cerebellum", cerebellum),
		("brainstem", brainstem),
		("hippocampus", hippocampus),
		("thalamus", thalamus),
		("amygdala", amygdala),
		("hypothalamus", hypothalamus),
		("gpt_net", gpt_net)
]
#
for model_name, model in models:
		model.sortModules()
#
def one_hot_encode(word, vocab_size, char_to_float):
		encoding = np.zeros((len(word), vocab_size), dtype=float)
		for i, char in enumerate(word):
				encoding[i, int(char_to_float[char])-7] = float(ord(char))
		return encoding

def one_hot_decode(encoding):
		decoded_word = ''
		for row in encoding:
				for val in row:
						if val != 0.0:
								decoded_word += chr(int(val))
		return decoded_word

vocab = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$%^&*()-=_+[]{}|;:\'",.<>?/\\`~ \n\t\v\r'
char_to_float = {char: float(i) for i, char in enumerate(vocab)}
vocab_size = len(vocab)

class BrainModel:
		def __init__(self, models):
				self.models = {name: model for name, model in models}
				self.shared_memory = np.zeros(100)
		def forward(self, sensory_input, previous_decision=None, social_input=None):
				thalamus_output = self.models['thalamus'].activate(sensory_input)
				emotional_state = self.models['amygdala'].activate(thalamus_output)
				autonomic_state = self.models['hypothalamus'].activate(thalamus_output)
				mirror_output = self.models['mirror_net'].activate(social_input if social_input is not None else np.zeros_like(sensory_input))
				if previous_decision is not None:
						self.shared_memory = self.models['hippocampus'].activate(previous_decision)
				decision_input = np.concatenate([thalamus_output, emotional_state, self.shared_memory, mirror_output, autonomic_state])
				decision = self.models['prefrontal_cortex'].activate(decision_input)
				executive_functions = self.models['frontal_lobe_net'].activate(decision)
				temporal_output = self.models['temporal_lobe'].activate(self.shared_memory)
				basic_functions = self.models['brainstem'].activate(thalamus_output)
				language_output = self.models['gpt_net'].activate(executive_functions)
				acc_output = self.models['acc_net'].activate(decision)
				language_input = np.concatenate([acc_output, decision])
				language_processing = self.models['brocas_area_net'].activate(language_input)
				comprehension = self.models['wernickes_area_net'].activate(language_processing)
				visual_input = self.models['visual_cortex'].activate(thalamus_output)
				final_input = np.concatenate([decision, executive_functions, temporal_output, basic_functions, language_output, acc_output, language_processing, comprehension, visual_input])
				final_decision = self.models['cerebral_cortex'].activate(final_input)
				return final_decision, comprehension
#
brain = BrainModel(models)
#
sensory_input = np.random.rand(400)
previous_decision = np.random.rand(250)
social_input = np.random.rand(500)
#
while True:
		user_text = input("Enter your input: ")
		one_hot_user_text = one_hot_encode(user_text, vocab_size, char_to_float)
		sensory_input = one_hot_user_text.flatten()
		social_input = np.concatenate([social_input, one_hot_user_text.flatten()])
		if len(sensory_input) < 400:
				sensory_input = np.pad(sensory_input, (0, 400 - len(sensory_input)), 'constant')
		elif len(sensory_input) > 400:
				sensory_input = sensory_input[:400]
		if len(social_input) < 500:
				social_input = np.pad(social_input, (0, 500 - len(social_input)), 'constant')
		elif len(social_input) > 500:
				social_input = social_input[:500]
		if previous_decision is None:
				previous_decision = np.zeros(250)
		final_decision, comprehension = brain.forward(sensory_input, previous_decision, social_input)
		print(final_decision, comprehension)
		previous_decision = np.pad(final_decision, (0, max(0, 250 - len(final_decision))), 'constant')[:250]

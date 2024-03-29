# What are these files for?

Templating AI of human heads. AI needs something to work with and off of. Since it doesn't have a biological corpus to work with, it needs something to emulate, simulate, and
 spoof til it matches human brain function per region of the brain. It is to have a model that dynamically reshapes and evolves as we do, but with some shortcuts.

# Why 3d neural network layers? What benefit does that give you?

It allows for 3 dimensions for thought. For instance, neurons pertaining to spatial, temporal, and other sensory information. Which, when it comes to AI development is essential for
 reasoning, imagining and otherwise a spatial dimension to work with like we use. Our cells have 3 coordinate points, so it is reasonable to use over excessive transformer models.
If you have say 12 transformers that use 2d grid hashing of information, you'll be their a lot longer than with 4 that have 3 dimensions to them. Meaning the 2d counterpart is
 limited to 2d maps of reality fed in bit by bit as 1d arrays of information.

For a human mind, 2-dimensions is not enough to describe how the world around them and how they imagine things. Because of that limitation, 3-dimensional layers would be more suited. As the 3-dimensions would be more realistic to come closer to human mental modeling. I then have to use it with brian2 & rdkit to make use of smiles and test for groups with different transmitters & receptors they're meant to imitate. 

However, I will give the docker container I would recommend building off of. 

https://hub.docker.com/repository/docker/de3343/ai_mods_py3.10/general

read the readme for the tag name for the latest build as it is not "latest" as listed on docker hub.

If you need more information, look to, https://chat.openai.com/share/1c475f8e-1fec-4348-9318-e65cdac45882, for frame of reference.

# Why does AI require brain emulation and data from a live host?

Because real data works better than faked or stolen data, token parsing isn't enough, and why not how to learn faster from a AI version that learned from you [you]? 

Among other things, AI should make use of eeg or the like data in live time with live, asynchronous training on the new data.
Now, resetting the dataset used for training while training it on something new allows for more degrees of freedom, keeping the window of error small.
In some cases, it will come down to manually tuning the knobs and dials to get EPICAC level of conscious machine then retooling the functions to automate this process a little without optimizing.

# Which `updated_main.py` should I use?

The default is fine to use as well as the latest one containing 3 layers each with their own 3-d orient. The ``dim3_neuronlayer.py`` has also been updated so it can optimally run
 based on the number of heads. This model could be adjusted for your needs. For instance, it could be configured to show what it is doing with the data as it attempts to respond.

# Why are you using the python ``pickle`` module?

This is to save the networks that are stable, with layers that work even for training the model by layer they're attached to.

The reason to use pickle is to save and reload the entire model. This may be temporary as the pybrain3 xml library for writing and reading xml snapshots is broken currently.
 Until a patch is used to fix it, the ``pickle`` module will have to do for now. 

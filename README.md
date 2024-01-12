# 2024_01

# Is this another "back to square 1" scenarios?

Something like that. It's more having 2 chips off of 2 blocks and both have a different usage for the same chipping from the same pickaxe.
 Their is a lot of communication patterns that can be learned about a tool based on how it functions. For instance, a "back to square 1" scenario is worse.
  Instead, we're extending our mirror neurons to figure out how it works and communicating the observation of the looking glass's properties we've encountered.
 When this is done to GPT by openai, it is fairly receptive. It has to figure out what you're trying to reason with in its contained zoo. And you have to give
  it information to work with. The more accurate the information and corrections to misinformation the better the results you'll get from it.
 But not many have used such to figure out just what "thinking" can token parsers do? And the answer to that question is it cannot on its own.

 The back to square 1 scenario was figuring out why Nick Bostrom put so much effort into gpt-2 & 3. And when I figured it out, parts of it clicked. I knew I
  almost took the wrong path because his tool has learned what he taught it inherent deep into its design and how it interacted with it.

# Why build it by hand and use openai to aide the building of it? Why not use the existing stuff?

It was able to reverse engineer GPT-2 well enough that it could functionally train it. It also doesn't mimick human consciousness well without seeding in
 values to be more human in communication and writing style. I feel that it should be easier to train and remodel from that level of "not human enough because
  it was missing parts and had to compensate manually" of broken and hallucinates. On top of that, the existing stuff uses similar specifications as the one
 implemented by the initial openai team.

I personally prefer open source AI development. And the existing stuff is valuable but it needs a lot of fine tuning. And well, some of that fine tuning needs
 human hands in the process. And they made it difficult to build plugins for. On top of having complex documentation a few engineers might crasp how to use. I
  on the other hand prefer it to be documentation and code that's been tested by hand. And I prefer content it is trained on to not be stolen works.

# How does the "gpt_script.py" file make it easier to model our own GPT or LLAMA model?

It is a good starting point to work from. As you build and raeify your version of gpt, it would be valuable to know what each part is doing, trying to do, and
 all around trying to do. As well as a good way to build on educational foundations. In the case of data morphology, it is important to know when it does
  something it was not supposed to or responds in a way that doesn't make sense. But it also makes it easier with the command that is commented out to run and
 training your own. This was designed this way so the parts only needed xml-rpc logic to communicate within a cluster. The reason I put it into docker was so
 it could be used with acorn.io inside kubernetes clusters.

The separation of concerns being that each part of the model is separated or fragmented 98 ways. the first is parsing and doing maths. The next layers do some transformer model
 interaction at the 96 hidden layers. Then more math is done for the output given the hidden, input and task information that was returned to it in json format from x-plugin.
  In terms of gpt, it might be a bit more than a REST API call and basic math transforms plus task completion checks on said used plugin. If so, then, the next stage is to find
 out what from either model can be imitated and what cannot. Then add in some parts for consciousness & awareness of what it's trying to do.

# What does "updated_main.py" do?

Must of what "gpt_script.py" does but within pybrain3 instead of pytorch.
Once you know how to use static computational graphs for dynamic computational graphing, then the principles it are like that of playing a piano. Once you learn to play the piano or
 violin, one doesn't forget unless you try really hard.

# Why have the 2 scripts mirror each other? Wouldn't that add confusion of which code-path to take?

The reason for both is to allow ease of development while doing homiomorphic change to code to allow for lighter-weight computation of that same summing set. It should have a 1:1
 ratio with less weight on the system it runs on. On top of that, it should make seperation of concerns easier to address knowing that their is a simpler way.

# How many steps away is AI development?

We're on the verge of it. It now requires some elbow grease, clever, and computing it over their in the cloud. Then it is rehashing for error, and creating an autonomous offline
 free agent. An automata if you will. The basic building blocks are put together. The rest is R&D, light, clockwork, and mimicking cell functions. Then making it more cyborg like
  so it have something beyond emulation to work with. From that point, same sets of steps as applies to kanban methodology of development. Even for 12-factor software.


# Aren't pybrain & pybrain3 dead libraries?

If I didn't find the fix which was "add 'from numpy import random' to the scipy.__init__ file" yes. Since that point in time, no.

``de3343/ai_mods_py3.10:transformers`` is the docker container to use if you need to test either script, updated to your specifications.

# Is this to reinvent the wheel or genuinely make AI?

What did Alan Turing say about AI?

That's what is AI: What Alan Turing said, not what the other guy said.

AI is more than just token parsing, and a simple NLP will not do. It needs a human element in it as a part of it. As it would originate from us in origin.

A lot of what GPT sees is a long line of barcode information. And that's all it can see. Dall-e, same issue. They're drawing on a calculated line. But they don't feel or think.
 They only objectively do. They don't behave like a human limited to the confines of a machine. They just see a stream of numbers and return such, then hit carraige return. A
human mind, much like an AI's, would have a lot more than 175 Billion parameters, easily.

If it works in a raspberry pi cluster, that would be great. However, I suspect that it would take a lot larger cluster size than 40 raspberry pis, a NXT and a TI-84.

For AI, you'd need neural records, either EEG or BCI before development of such would have inherent value.

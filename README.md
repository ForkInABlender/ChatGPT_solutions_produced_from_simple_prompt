# Why minecraft & Python3.10?

Minecraft allows for sandboxed interactions. Python makes many points in execution of code very flecible as a c++ based runtime. 

Jython and the like weren't good enough for the type of inter-language inter-runtime perefusion of execution. Jpype had a better interface. Minecraft has many silly features.

But it has preprogrammed "decisions" for its npc's. The reason I don't condone AI based video game development is video game development underminds what is needed for a human mind 
 to express itself creatively. Minecrafts' initial development team didn't limit all of what the user can do. It limited what couldn't be done. 

What it would allow for AI development as a video game instead is a space for a human mind to make use. Otherwise, the human mind in AI terms would be floating in free space.
 unable to do anything about it while aware of itself and what it learned per say. While minecraft is a video game, I would not build AI into it. But around and making use of it,
 sure. Make use of minecraft with AI doing UI feedback looping over the context it has to operate with, different story. As that is making use of minecraft is meant to be temporary.

The goal is to use it to aide in modeling neural networks that can accept a 3-d environment and have similar periferal vision as humans do while playing the game. Later on, this can
 be expanded to having a network modeled for "imagining" what next it can do. Then give it something more to work with.


# Why GTK and OpenGL in addition in python?

Ease of development. The only thing that is being done different is how it is made use of. OpenGL in python allows me to hide away the rendering that the AI would do on a GPU.
 For instance, you're using an intel gpu on a dell office computer. I personally don't envy hopping through code bases and drivers for a simple answer stairing at me.
 PyOpenGL worked for using OpenGL targeting a mesa driver, where the gtk framework would work on other builds. 

# What is GridShift Encryption for?

It is in this repo to make use of it in AI. The point is to train it on encoded responses and return similar. Then decode and return response. The reason I find it important for AI
 is to capture generated and pretrained responses. The end goal is to make use of needed tokens rather than just make use of large numbers & large numerical sets. Instead, small
 changes are required for better sets of responses. So, part of my thinking was encrypt & tokenize then train and detokenize and decode the response. The gridshift encryption 
 can also be applied to AI oriented neural network modeling for dynamic reshape, assuming the layer sizes match or follow a custom grid-wrap rule. This would also allow for perturbing your network somewhat to your desired outcome. In theory, it would allow for also shifting neurons around from any layer so long as you knew when to shift the neurons or where to keep responses encrypted based on the input and shift number & ascii tokens used. In the available example, it does a 10 by 10 of common to use characters for English speakers. This later can be adapted to other languages as well, but, it should provide a basic outline to work with.

# What does Acorn.io have to do with building an AI?

Acorn.io allows a lot of flexibility with how you structure, constrain, control, or operate kubernetes. Even from their source source code. What it gives python developers is the 
 option to use something like tiktoken, pybrain, luigi, and pynisher with flask. Acorn also allows for ease of definition of what and how it should conditionally respond.

The point of using acorn.io is mostly to cut the time needed to train each copy of the same container image of the pybrain3/tiktoken version. So that development time doesn't rely 
 on hardware like gpu heavy compute. Instead, make it simple, lightweight, and able to run anywhere. In my case, it isn't about putting an AI into a video game, but giving it a 
 framework to run around in so the human it is emulating doesn't go insane. The point of it being able to make use of minecraft is to have something that is psuedo realistic
 enough so it can interact with A environment which gives it freedom of mobility and communicationand no limitations on thinking. Instead of trapping a human mind, the goal is to
 template from this point til it matches that of the human brain with no goals for it to be limited to minecraft game world data. ( A world seed if you will ... )

# Why not do it proprietarily? Isn't open source unsafe?

Proprietary development is not the point. The reason it is being done in this manner as open source is that proprietary development focused and hyperfocused on token parsing, which isn't AI so much as poorly designed autocorrect with the intent to fool. I do not, for that reason, support OpenAI's lack of critical thinking on AI development. Even Google doesn't
 have the chops for AI development in the real world. Same is true for Microsoft spending billions in infra that is only impeding development. However, Google is a lot closer than many have come. While google isn't proprietary inherently, they do follow similar bending of the knee this way problematically. 

Open source development, even from what pybrain/pybrain3's initial developers created, proved to be more functional than modern day proprietarily built libraries and software. What open source development of AI allows for is the same scope analysis to the problem set NASA uses when finding the root causes of a problem.


While still unsafe, tested and mildly unsafe is better than absolute panick, which is what openai has as a code base. At best it is a good template to work from, but, lacks human qualities that would be considered inherent by Alan Turings' definition of AI, much like Isaac Asimovs'.

# What would happen if someone used this code to make an AI of you?

They'd get trolled by the model once they'd have succeeded. Imitation is not flattery in that case, and the results they get are what they get. Now, I will not list the exceptions to that rule because that is for me and that person the exception applies to  to know. And if they use it to lie, then the results will be just that based on the input. GIGO law applies to it and why you did that with it. And that's if you force it you play your lie out as truth.

# Aren't you going to give a how to guide with the readme?

Why? It's your head. Organize your head how you wish. Their is no one way on how to think. Either you know how to use it, or you do not. And you'd know how to use it if you read the source material it is built on. For instance, pybrain3 based network layers, networks, etcetera, are mostly for static computational graphs but can still do dynamic computational graph representation similar to what gpt-4 uses pytorch for under the hood. tiktoken is good for token parsing. 

AI development has several purposes:
 * problem solve
 * imitate a human mind represented as dynamically reshaping networks
 * Do more than token parse & actually comprehend beyond basic training
 * imitating human consciousness via rdkit+brian2+pybrain3 based neural network ( not for the purposes of embedding into video games as that violates the individual it is imitatings' copyright rights including their human rights )

The reasoning for such is that it must also use real values for biological functions, sensory inputs, and chemical responses a human does. A video game environment would only drive the person that AI is imitating insane, as they could never functionally escape the zoo you put them in. So, it must be able to see where and what it is operating and outside of itself, a reprinted copy of its body from genetics gotten from its conscious mind. While this is a very noisy process, due to the static, resamplings of data will be required to do noise cancelation. The point of the noise cancelation will be to prevent random signals from eeg data or bci data from causing misalignment and prevent one from missing anything specific for a 1:1 match of their mind from the raw data.

# Are you planning on making snapshots of each using xml documents?

Yes. The key to valid and savable layers is the ``self.name`` has to be manually assigned in any custom layer, connection, etcetera, that you design from here. As long as it has
 a name to be saved by. 

# why `udocker` and why Android OS? Why not use the apps on the android market and just use ssh?

It makes it easier and saves me ( and others ) development time. Plus ssh is a pain in the ass, and the chiefs that make use of it, more power too them. I, however, found ssh to 
be a stumbling block impeding development. It also had to be run offline before it would be useful in many clusters beyond a cell phone. On top of that, it saved me time and development cost. In addition, it was easier to consider the minimum build that would work for systems using proot-distro. ( I can't speak for IOS and apple software development, however, they may have similar enough binaries that it would also work and run on an iphone. )

``alias docker="udocker"`` may come in handy in the proot-distro of an ubuntu install.

While limited, it can run docker containers within itself with a docker-in-docker solution. Meaning it can still be useful in other ways while passing the docker socket to the udocker instance of that [D-in-D] solution to the limitations of udocker.

Sadly, I may have made it every OS compatible......
Forunately, this is open source code. So, "yay, it can run everywhere including a cell-phone". Unfortunately, that was found on a very high up on a steep gradient. So it is a very (vuuury) hard thing to build.


But, that being said, the code within this repo branch should make it easier. if you're having trouble getting your copy of the docker container working correctly with the code:

``docker run -it --rm -v $PWD/:/app/ de3343/ai_mods_py3.10:transformers sh -c "mv /string.py /usr/local/lib/python3.10/string.py; python3 /app/updated_main.py"`` within the cloned repo on your machine or phone. Because the code within it is python3.10 compatible, anything kept compatible in 3.10+ will still work. If you need help getting it to work on your android phone, ask chatgpt and check the docs for your android version. If your android OS is a x86_64 build, the default docker build should work.

# Why would you use GPT to aide in the modeling of an AI? Why not just use the GPT model and not care?

A human mind does more than just parse tokens. Or as Terry Pratchett put it in his discworld series "from the outside, infinity is blue; it is only black while you're standing inside it". And a human mind as dynamic neurochemical responses to its environment. It equally needs to emulate or imitate, if not be part of, a living host. GPT can also be gainsayed in a way that a human can't be coerced. 

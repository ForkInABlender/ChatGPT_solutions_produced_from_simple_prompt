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

# Why pybrain3?

It kept the development of a framework for an AI to be built off of simple. It also made it easier to integrate bits around it that an AI would use when housed inside its toolset.
 Amongst other things, it made it easier due to its initial developer's work on the pybrain3 neural network programming ecosystem. Torch, Tensorflow, and Transformer modeling
  architecture optimized away the needed fine grain values that would apply to a sentient being thinking about the inputs it would be about to or currently responding to.

Excessive optimization produces deadwood and adds bloatware that would have initially killed the project.

The projects' aim is to imitate the human head based on language input, visual input, and auditorial input, plus some abstractions needed for neurochemical emulation as close to
 the dynamic a human head uses. For reference look to the folder "pybrain networks" and "MAML". Under that folder, you'll find bits needed to construct an AI before needing EEG, BCI,
  fMRI (non-intrusive, small scale DIY & the like) and extrapolate backward from that data.

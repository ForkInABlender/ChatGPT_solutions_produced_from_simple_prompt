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
 have the chops for AI development in the real world. Same is true for Microsoft spending billions in infra that is only impeding development. While google isn't proprietary inherently, they do follow similar bending of the knee this way problematically. 

Open source development, even from what pybrain/pybrain3's initial developers created, proved to be more functional than modern day proprietarily built libraries and software. What open source development of AI allows for is the same scope analysis to the problem set NASA uses when finding the root causes of a problem.

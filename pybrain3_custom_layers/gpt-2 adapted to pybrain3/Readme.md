# GPT-2 now uses numpy? What gives?

The version here uses numpy & Scipy instead of tensorflow. The choice was for simplified remodeling later so it had the proper layers, layout, and basic functionality.

``pybrain3`` allows for similar imitation but further reduces the potential for error. Beyond that, test as you'd like as this version was not written by a proprietary company.

As the model I am building, while it touches on similar approach, it requires the neurochemical and biochemical simulator bits and eeg/bci data for genuine AI to be modeled.
Token parsing is not good enough as token parsing, as gpt-2, 3, & 4, as presently built by openAI developers, is nothing more than a suffisticated form of the Eliza project, designed to defeat the purpose of the turing test with the intent on fooling the user with clever symbolic references. 

Due to openAI's, and many others who failed to pay attention to what Alan Turing said, Alan Turing's version of AI, AI will be built by open source developers like myself.

# why does that differential matter? Isn't attention all that is needed?

Sure, if one is developing it, but, that lacks emotional reasoning and argues that inferencing and symbolic references is all that's needed. One question: Where is it's brain in
 the process? No where, so it can't be called "sentient" or "artificial intelligence". By that notion, you might as well call your coffee table or smart fridge "AI". Or for that matter, your dresser. Is your dresser or coffee table sentient and harboring artificial intelligence? I'm pretty sure you'd agree that the last 2 questions are an indication of why assuming token parsing and attention is not all that is needed for AI.

# So the hype about "AI" and what's on the market currently isn't AI?

Correct. It's as useful as the Eliza projects in conversation and "reasoning" \*\*ahem\*\*. 

# How long before your model is an AI based on your pybrain3 code?

Depends, really. I am aiming for a AI development job. So, It will take significantly longer than initially anticipated. It also depends on the training it gets, data used, human interactive testing, etc. And when new layers get added it isn't simply transfering the data and recalibrating the old and new layer. It is more or less retraining the entire model.

I'd say it would take someone from this point on 3 months flat with a kubernetes cluster at their disposal. If they're using k3d/l3s on their laptop or cell phone, probably a bit longer. Even then, it has to go through the same tests a normal human would ( think of how a neurologist & psychologist break down problems ). It would take just as long to to train. For instance, it needs to be able to emulate a GABA, CB1, & CB2 junctors most if not all living beings with a brain would. In this case, I'm aiming for it to have enough
 to emulate/simulate til it isn't a simulated or emulated version and proves to have a 1:1 ratio to human nerve cells work. 

# Then what, have it map their neural patterns to genetic templates?

Yes, actually. It's one function of the brain. Our genetics in our heads change and why should the AI templates do any different in that respect? Sure, there are bound to be errors, but those errors will be easier to pinpoint as human error rather than assuming it to be natural phenomenon. If a person changes IQ ranges, the genes governing that reshape due to the alternative which would be death (in this case, it decreases rather than increases..."death"....).

However, this would also give more insight to how diseases, disorders, and the like. It would be more useful to know with a higher degree of precision on what is failing and where to look rather than guessing.

By the same token, AI based off of a human mind including the neurochemical and biological functions at play for it would be just as human. In realistic terms, it should have a genetic map using the same 4 base pairs we use.

# How would you know if you got it right?

There is only one way to find out :)
minimum going rate for genuine AI development by a developer such as myself: $1000/hr, minimum. Otherwise, the results you get are likely guess and check, and probably more agonizing than finding the one qualified with a fair going rate. 

 # Why does this work with ``de3343/ai_mods_py3.10:numba_with_fire_included``? Why can't I use torch, torchaudio, or torchvision like I used to in the prior image?

That was by design. Simplicity and flexibility. Torch is valuable but it is heavier, oversimplifies via optimization, versus simply not reinventing the wheel...

Instead use python3.10's native RPC libraries with the correct port and ip exposes so it can be able to use those functions between containers. While RPC can be an insecure way to make use of those libraries defined in the prior docker image, it is recommended over excessive optimization producing dead wood. Torch, in this case, is the dead wood. And I remove the dead would by throwing it into the code graveyard when I produce such deadwood code I believe could be repurposed. 

( Maybe like ``pybrain3`` it isn't fully dead, but, if you don't poke it with a stick, you won't be sure if it is deadwood or still a viable codebase to work with; and  ``pybrain3`` is 2008 old and still able to outperform current libraries by its developer's design choice. I only checked one thing about python3.6 string.py and python3.10 then hit run. Which brought ``pybrain3`` back to its feet and regrow its skin for analogies sake. )

While ``pybrain3`` is 12 years out of date, it still can do the dynamic computation and transformer logic wtih some clever/intelligent application of elbow grease. Even though its initial developers probably had not thought it possible to use a static computational graph to perform similar to a dynamic one. None the less, I figured out how it could be done in a practical sense.

# Why not use the modern day ML libraries for python over ``pybrain/pybrain3``? Isn't the point efficiency?

Yes, it is about efficiency, but, no, modern day libraries use a lot of boilerplate (dead wood) that is poorly maintained which I can not accept usage of. Poor design by abstracting away the complex parts, making using the modern-day ML libraries inefficient and bloatware I typically discard.

# what is a code graveyard?

Often a code grave yard is a graveyard full of tomb stones with QR codes on them. It is also often referred to as a place where code goes and is expected to have no use or little to no potential for refactor due to wrong approach to a problem ( like a trash bin for chunking notes you'd crumple up and toss rather than burn because maybe you missed something; but maybe one thought it the dumb way and simply overlooked a small detail said prior note or "code" did not account for; that'd be reason to turn the trash bin upside down knowing it contained only crumpled up notes )

This is an example of my [graveyard](https://github.com/ForkInABlender/Project_Graveyard/) for instance. I rarely make releases of the deadwood accumulated ( or crumpled up notes I tossed into the rubbish bin )

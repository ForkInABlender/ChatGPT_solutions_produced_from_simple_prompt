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

# what is a code graveyard and why should such be imitated as you've laid yours out?

Often a code grave yard is a graveyard full of tomb stones with QR codes on them. It is also often referred to as a place where code goes and is expected to have no use or little to no potential for refactor due to wrong approach to a problem ( like a trash bin for chunking notes you'd crumple up and toss rather than burn because maybe you missed something; but maybe one thought it the dumb way and simply overlooked a small detail said prior note or "code" did not account for; that'd be reason to turn the trash bin upside down knowing it contained only crumpled up notes )

This is an example of my [graveyard](https://github.com/ForkInABlender/Project_Graveyard/) for instance. I rarely make releases of the deadwood accumulated ( or crumpled up notes I tossed into the rubbish bin )

The reason to imitate the process is it allows for deadwood to be looked at in a batch mannerism for cross inferencing against. It saves the failed attempts and one can give details of what resulted in some of it, or it is meant to stay truly dead. If someone does stop the part where the developer spun their wheels and gave up, it would be good for code review even by those not on any tech team you work with as well as those you do. What it also allows the open door for is basically seeing where one did the xy-problem analysis wrong to the part they're aiming to solve for. If you can't spot what is done wrong or why the calculation and the result aren't adding up, you might not be able to find it & end up giving up entirely. Where putting it into similar format of a code-graveyard, isn't giving up so much as passing the baton. 


# Are their any vulnerabilities on the container ``#``?

Yes.


<pre><b>114 vulnerabilities found in 43 packages</b>
  <font color="#FFFFFF">UNSPECIFIED</font>  <font color="#FFFFFF">2</font>   
  <span style="background-color:#FFD787"><font color="#000000">LOW</font></span>          <font color="#FFD787">79</font>  
  <span style="background-color:#FFAF5F"><font color="#000000">MEDIUM</font></span>       <font color="#FFAF5F">22</font>  
  <span style="background-color:#D75F5F"><font color="#000000">HIGH</font></span>         <font color="#D75F5F">11</font>  
  CRITICAL     0  </pre>
<pre><font color="#FF8787"><b>pkg:</b></font>deb/debian/libde265@1.0.11-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/postgresql-15@15.3-0+deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/gnutls28@3.7.9-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/expat@2.5.0-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>pypi/pillow@10.1.0
<font color="#FF8787"><b>pkg:</b></font>deb/debian/imagemagick@8:6.9.11.60+dfsg-1.6?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/tiff@4.5.0-6?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/libwmf@0.2.12-5.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/openssh@1:9.2p1-2+deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/curl@7.88.1-10+deb12u4?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/libgcrypt20@1.10.1-3?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/systemd@252.17-1~deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/mariadb@1:10.11.4-1~deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/bluez@5.66-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/dav1d@1.0.0-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>pypi/jinja2@3.1.2
<font color="#FF8787"><b>pkg:</b></font>deb/debian/openjpeg2@2.5.0-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/binutils@2.40-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/patch@2.7.6-7?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/openldap@2.5.13+dfsg-5?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/glibc@2.36-9+deb12u3?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/git@1:2.39.2-1.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/shadow@1:4.13+dfsg1-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/openssl@3.0.11-1~deb12u2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/m4@1.4.19-3?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/tar@1.34+dfsg-1.2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/unzip@6.0-28?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/elfutils@0.188-2.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/gcc-12@12.2.0-14?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/apt@2.6.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/nghttp2@1.52.0-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/jbigkit@2.1-6.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/openexr@3.1.5-5?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/libpng1.6@1.6.39-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/krb5@1.20.1-2+deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/pixman@0.42.2-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/libxslt@1.1.35-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/util-linux@2.38.1-5+b1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/util-linux@2.38.1-5?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/gnupg2@2.2.40-1.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/coreutils@9.1-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/libheif@1.15.1-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12
<font color="#FF8787"><b>pkg:</b></font>deb/debian/perl@5.36.0-7?os_distro=bookworm&amp;os_name=debian&amp;os_version=12</pre>
<pre>
<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#D75F5F"><font color="#000000">   5H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   2M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>libde265 1.0.11-1</b>
pkg:deb/debian/libde265@1.0.11-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#D75F5F">✗ HIGH</font> CVE-2023-49468
      https://scout.docker.com/v/CVE-2023-49468
      Affected range : &lt;1.0.11-1+deb12u2  
      Fixed version  : 1.0.11-1+deb12u2   
    
    <font color="#D75F5F">✗ HIGH</font> CVE-2023-49467
      https://scout.docker.com/v/CVE-2023-49467
      Affected range : &lt;1.0.11-1+deb12u2  
      Fixed version  : 1.0.11-1+deb12u2   
    
    <font color="#D75F5F">✗ HIGH</font> CVE-2023-49465
      https://scout.docker.com/v/CVE-2023-49465
      Affected range : &lt;1.0.11-1+deb12u2  
      Fixed version  : 1.0.11-1+deb12u2   
    
    <font color="#D75F5F">✗ HIGH</font> CVE-2023-27103
      https://scout.docker.com/v/CVE-2023-27103
      Affected range : &lt;1.0.11-1+deb12u1  
      Fixed version  : 1.0.11-1+deb12u1   
    
    <font color="#D75F5F">✗ HIGH</font> CVE-2023-43887
      https://scout.docker.com/v/CVE-2023-43887
      Affected range : &lt;1.0.11-1+deb12u1  
      Fixed version  : 1.0.11-1+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-47471
      https://scout.docker.com/v/CVE-2023-47471
      Affected range : &lt;1.0.11-1+deb12u1  
      Fixed version  : 1.0.11-1+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-27102
      https://scout.docker.com/v/CVE-2023-27102
      Affected range : &lt;1.0.11-1+deb12u1  
      Fixed version  : 1.0.11-1+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#D75F5F"><font color="#000000">   3H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   2M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>postgresql-15 15.3-0+deb12u1</b>
pkg:deb/debian/postgresql-15@15.3-0+deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#D75F5F">✗ HIGH</font> CVE-2023-5869
      https://scout.docker.com/v/CVE-2023-5869
      Affected range : &lt;15.5-0+deb12u1  
      Fixed version  : 15.5-0+deb12u1   
    
    <font color="#D75F5F">✗ HIGH</font> CVE-2023-39417
      https://scout.docker.com/v/CVE-2023-39417
      Affected range : &lt;15.5-0+deb12u1  
      Fixed version  : 15.5-0+deb12u1   
    
    <font color="#D75F5F">✗ HIGH</font> CVE-2024-0985
      https://scout.docker.com/v/CVE-2024-0985
      Affected range : &lt;15.6-0+deb12u1  
      Fixed version  : 15.6-0+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-5870
      https://scout.docker.com/v/CVE-2023-5870
      Affected range : &lt;15.5-0+deb12u1  
      Fixed version  : 15.5-0+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-5868
      https://scout.docker.com/v/CVE-2023-5868
      Affected range : &lt;15.5-0+deb12u1  
      Fixed version  : 15.5-0+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-39418
      https://scout.docker.com/v/CVE-2023-39418
      Affected range : &lt;15.5-0+deb12u1  
      Fixed version  : 15.5-0+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#D75F5F"><font color="#000000">   1H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>gnutls28 3.7.9-2</b>
pkg:deb/debian/gnutls28@3.7.9-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#D75F5F">✗ HIGH</font> CVE-2024-0553
      https://scout.docker.com/v/CVE-2024-0553
      Affected range : &lt;3.7.9-2+deb12u2  
      Fixed version  : 3.7.9-2+deb12u2   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-5981
      https://scout.docker.com/v/CVE-2023-5981
      Affected range : &lt;3.7.9-2+deb12u1  
      Fixed version  : 3.7.9-2+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2024-0567
      https://scout.docker.com/v/CVE-2024-0567
      Affected range : &lt;3.7.9-2+deb12u2  
      Fixed version  : 3.7.9-2+deb12u2   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#D75F5F"><font color="#000000">   1H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <font color="#FFFFFF">   1? </font> <b>expat 2.5.0-1</b>
pkg:deb/debian/expat@2.5.0-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#D75F5F">✗ HIGH</font> CVE-2023-52425
      https://scout.docker.com/v/CVE-2023-52425
      Affected range : &gt;=2.5.0-1  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-52426
      https://scout.docker.com/v/CVE-2023-52426
      Affected range : &gt;=2.5.0-1  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFFFFF">✗ UNSPECIFIED</font> CVE-2024-28757
      https://scout.docker.com/v/CVE-2024-28757
      Affected range : &gt;=2.5.0-1  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#D75F5F"><font color="#000000">   1H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>pillow 10.1.0</b>
pkg:pypi/pillow@10.1.0

    <font color="#D75F5F">✗ HIGH</font> CVE-2023-50447 [Improper Control of Generation of Code (&apos;Code Injection&apos;)]
      https://scout.docker.com/v/CVE-2023-50447
      Affected range : &lt;10.2.0                                       
      Fixed version  : 10.2.0                                        
      CVSS Score     : 8.1                                           
      CVSS Vector    : CVSS:3.1/AV:N/AC:H/PR:N/UI:N/S:U/C:H/I:H/A:H  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   3M </font></span> <span style="background-color:#FFD787"><font color="#000000">  13L </font></span> <b>imagemagick 8:6.9.11.60+dfsg-1.6</b>
pkg:deb/debian/imagemagick@8:6.9.11.60+dfsg-1.6?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-5341
      https://scout.docker.com/v/CVE-2023-5341
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-34151
      https://scout.docker.com/v/CVE-2023-34151
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-1289
      https://scout.docker.com/v/CVE-2023-1289
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-3428
      https://scout.docker.com/v/CVE-2023-3428
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-34152
      https://scout.docker.com/v/CVE-2023-34152
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-1906
      https://scout.docker.com/v/CVE-2023-1906
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2022-1115
      https://scout.docker.com/v/CVE-2022-1115
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2021-3610
      https://scout.docker.com/v/CVE-2021-3610
      Affected range : &lt;8:6.9.11.60+dfsg-1.6+deb12u1  
      Fixed version  : 8:6.9.11.60+dfsg-1.6+deb12u1   
    
    <font color="#FFD787">✗ LOW</font> CVE-2021-20311
      https://scout.docker.com/v/CVE-2021-20311
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-15607
      https://scout.docker.com/v/CVE-2018-15607
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-7275
      https://scout.docker.com/v/CVE-2017-7275
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-11755
      https://scout.docker.com/v/CVE-2017-11755
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-11754
      https://scout.docker.com/v/CVE-2017-11754
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-8678
      https://scout.docker.com/v/CVE-2016-8678
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2008-3134
      https://scout.docker.com/v/CVE-2008-3134
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    
    <font color="#FFD787">✗ LOW</font> CVE-2005-0406
      https://scout.docker.com/v/CVE-2005-0406
      Affected range : &gt;=8:6.9.11.60+dfsg-1.6  
      Fixed version  : <font color="#FF8787">not fixed</font>               
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   3M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>tiff 4.5.0-6</b>
pkg:deb/debian/tiff@4.5.0-6?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-41175
      https://scout.docker.com/v/CVE-2023-41175
      Affected range : &lt;4.5.0-6+deb12u1  
      Fixed version  : 4.5.0-6+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-40745
      https://scout.docker.com/v/CVE-2023-40745
      Affected range : &lt;4.5.0-6+deb12u1  
      Fixed version  : 4.5.0-6+deb12u1   
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-3576
      https://scout.docker.com/v/CVE-2023-3576
      Affected range : &lt;4.5.0-6+deb12u1  
      Fixed version  : 4.5.0-6+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   2M </font></span> <span style="background-color:#FFD787"><font color="#000000">   2L </font></span> <b>libwmf 0.2.12-5.1</b>
pkg:deb/debian/libwmf@0.2.12-5.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2009-3546
      https://scout.docker.com/v/CVE-2009-3546
      Affected range : &gt;=0.2.12-5.1  
      Fixed version  : <font color="#FF8787">not fixed</font>     
    
    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2007-3996
      https://scout.docker.com/v/CVE-2007-3996
      Affected range : &gt;=0.2.12-5.1  
      Fixed version  : <font color="#FF8787">not fixed</font>     
    
    <font color="#FFD787">✗ LOW</font> CVE-2007-3477
      https://scout.docker.com/v/CVE-2007-3477
      Affected range : &gt;=0.2.12-5.1  
      Fixed version  : <font color="#FF8787">not fixed</font>     
    
    <font color="#FFD787">✗ LOW</font> CVE-2007-3476
      https://scout.docker.com/v/CVE-2007-3476
      Affected range : &gt;=0.2.12-5.1  
      Fixed version  : <font color="#FF8787">not fixed</font>     
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#FFD787"><font color="#000000">   3L </font></span> <b>openssh 1:9.2p1-2+deb12u1</b>
pkg:deb/debian/openssh@1:9.2p1-2+deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-51385
      https://scout.docker.com/v/CVE-2023-51385
      Affected range : &lt;1:9.2p1-2+deb12u2  
      Fixed version  : 1:9.2p1-2+deb12u2   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-51384
      https://scout.docker.com/v/CVE-2023-51384
      Affected range : &lt;1:9.2p1-2+deb12u2  
      Fixed version  : 1:9.2p1-2+deb12u2   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-48795
      https://scout.docker.com/v/CVE-2023-48795
      Affected range : &lt;1:9.2p1-2+deb12u2  
      Fixed version  : 1:9.2p1-2+deb12u2   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-28531
      https://scout.docker.com/v/CVE-2023-28531
      Affected range : &lt;1:9.2p1-2+deb12u2  
      Fixed version  : 1:9.2p1-2+deb12u2   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>curl 7.88.1-10+deb12u4</b>
pkg:deb/debian/curl@7.88.1-10+deb12u4?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-46218
      https://scout.docker.com/v/CVE-2023-46218
      Affected range : &lt;7.88.1-10+deb12u5  
      Fixed version  : 7.88.1-10+deb12u5   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-46219
      https://scout.docker.com/v/CVE-2023-46219
      Affected range : &lt;7.88.1-10+deb12u5  
      Fixed version  : 7.88.1-10+deb12u5   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>libgcrypt20 1.10.1-3</b>
pkg:deb/debian/libgcrypt20@1.10.1-3?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2024-2236
      https://scout.docker.com/v/CVE-2024-2236
      Affected range : &gt;=1.10.1-3  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-6829
      https://scout.docker.com/v/CVE-2018-6829
      Affected range : &gt;=1.10.1-3  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>jinja2 3.1.2</b>
pkg:pypi/jinja2@3.1.2

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2024-22195 [Improper Neutralization of Input During Web Page Generation (&apos;Cross-site Scripting&apos;)]
      https://scout.docker.com/v/CVE-2024-22195
      Affected range : &lt;3.1.3                                        
      Fixed version  : 3.1.3                                         
      CVSS Score     : 5.4                                           
      CVSS Vector    : CVSS:3.1/AV:N/AC:L/PR:N/UI:R/S:U/C:L/I:L/A:N  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>systemd 252.17-1~deb12u1</b>
pkg:deb/debian/systemd@252.17-1~deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-7008
      https://scout.docker.com/v/CVE-2023-7008
      Affected range : &lt;252.21-1~deb12u1  
      Fixed version  : 252.21-1~deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>mariadb 1:10.11.4-1~deb12u1</b>
pkg:deb/debian/mariadb@1:10.11.4-1~deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-22084
      https://scout.docker.com/v/CVE-2023-22084
      Affected range : &lt;1:10.11.6-0+deb12u1  
      Fixed version  : 1:10.11.6-0+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>dav1d 1.0.0-2</b>
pkg:deb/debian/dav1d@1.0.0-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2024-1580
      https://scout.docker.com/v/CVE-2024-1580
      Affected range : &gt;=1.0.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#FFAF5F"><font color="#000000">   1M </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0L </font></span> <b>bluez 5.66-1</b>
pkg:deb/debian/bluez@5.66-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFAF5F">✗ MEDIUM</font> CVE-2023-45866
      https://scout.docker.com/v/CVE-2023-45866
      Affected range : &lt;5.66-1+deb12u1  
      Fixed version  : 5.66-1+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">  13L </font></span> <b>openjpeg2 2.5.0-2</b>
pkg:deb/debian/openjpeg2@2.5.0-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2018-20846
      https://scout.docker.com/v/CVE-2018-20846
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-16376
      https://scout.docker.com/v/CVE-2018-16376
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-16375
      https://scout.docker.com/v/CVE-2018-16375
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-17479
      https://scout.docker.com/v/CVE-2017-17479
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9581
      https://scout.docker.com/v/CVE-2016-9581
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9580
      https://scout.docker.com/v/CVE-2016-9580
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9117
      https://scout.docker.com/v/CVE-2016-9117
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9116
      https://scout.docker.com/v/CVE-2016-9116
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9115
      https://scout.docker.com/v/CVE-2016-9115
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9114
      https://scout.docker.com/v/CVE-2016-9114
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-9113
      https://scout.docker.com/v/CVE-2016-9113
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-10506
      https://scout.docker.com/v/CVE-2016-10506
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2016-10505
      https://scout.docker.com/v/CVE-2016-10505
      Affected range : &gt;=2.5.0-2  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   7L </font></span> <b>binutils 2.40-2</b>
pkg:deb/debian/binutils@2.40-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2023-1972
      https://scout.docker.com/v/CVE-2023-1972
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2021-32256
      https://scout.docker.com/v/CVE-2021-32256
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-9996
      https://scout.docker.com/v/CVE-2018-9996
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-20712
      https://scout.docker.com/v/CVE-2018-20712
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-20673
      https://scout.docker.com/v/CVE-2018-20673
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-18483
      https://scout.docker.com/v/CVE-2018-18483
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-13716
      https://scout.docker.com/v/CVE-2017-13716
      Affected range : &gt;=2.40-2   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   4L </font></span> <b>patch 2.7.6-7</b>
pkg:deb/debian/patch@2.7.6-7?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2021-45261
      https://scout.docker.com/v/CVE-2021-45261
      Affected range : &gt;=2.7.6-7  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-6952
      https://scout.docker.com/v/CVE-2018-6952
      Affected range : &gt;=2.7.6-7  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-6951
      https://scout.docker.com/v/CVE-2018-6951
      Affected range : &gt;=2.7.6-7  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    
    <font color="#FFD787">✗ LOW</font> CVE-2010-4651
      https://scout.docker.com/v/CVE-2010-4651
      Affected range : &gt;=2.7.6-7  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   4L </font></span> <b>openldap 2.5.13+dfsg-5</b>
pkg:deb/debian/openldap@2.5.13+dfsg-5?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2020-15719
      https://scout.docker.com/v/CVE-2020-15719
      Affected range : &gt;=2.5.13+dfsg-5  
      Fixed version  : <font color="#FF8787">not fixed</font>        
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-17740
      https://scout.docker.com/v/CVE-2017-17740
      Affected range : &gt;=2.5.13+dfsg-5  
      Fixed version  : <font color="#FF8787">not fixed</font>        
    
    <font color="#FFD787">✗ LOW</font> CVE-2017-14159
      https://scout.docker.com/v/CVE-2017-14159
      Affected range : &gt;=2.5.13+dfsg-5  
      Fixed version  : <font color="#FF8787">not fixed</font>        
    
    <font color="#FFD787">✗ LOW</font> CVE-2015-3276
      https://scout.docker.com/v/CVE-2015-3276
      Affected range : &gt;=2.5.13+dfsg-5  
      Fixed version  : <font color="#FF8787">not fixed</font>        
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   3L </font></span> <b>glibc 2.36-9+deb12u3</b>
pkg:deb/debian/glibc@2.36-9+deb12u3?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2023-6780
      https://scout.docker.com/v/CVE-2023-6780
      Affected range : &lt;2.36-9+deb12u4  
      Fixed version  : 2.36-9+deb12u4   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-6779
      https://scout.docker.com/v/CVE-2023-6779
      Affected range : &lt;2.36-9+deb12u4  
      Fixed version  : 2.36-9+deb12u4   
    
    <font color="#FFD787">✗ LOW</font> CVE-2023-6246
      https://scout.docker.com/v/CVE-2023-6246
      Affected range : &lt;2.36-9+deb12u4  
      Fixed version  : 2.36-9+deb12u4   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   2L </font></span> <b>git 1:2.39.2-1.1</b>
pkg:deb/debian/git@1:2.39.2-1.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2022-24975
      https://scout.docker.com/v/CVE-2022-24975
      Affected range : &gt;=1:2.39.2-1.1  
      Fixed version  : <font color="#FF8787">not fixed</font>       
    
    <font color="#FFD787">✗ LOW</font> CVE-2018-1000021
      https://scout.docker.com/v/CVE-2018-1000021
      Affected range : &gt;=1:2.39.2-1.1  
      Fixed version  : <font color="#FF8787">not fixed</font>       
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   2L </font></span> <b>m4 1.4.19-3</b>
pkg:deb/debian/m4@1.4.19-3?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2008-1688
      https://scout.docker.com/v/CVE-2008-1688
      Affected range : &gt;=1.4.19-3  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    
    <font color="#FFD787">✗ LOW</font> CVE-2008-1687
      https://scout.docker.com/v/CVE-2008-1687
      Affected range : &gt;=1.4.19-3  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   2L </font></span> <b>shadow 1:4.13+dfsg1-1</b>
pkg:deb/debian/shadow@1:4.13+dfsg1-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2019-19882
      https://scout.docker.com/v/CVE-2019-19882
      Affected range : &gt;=1:4.13+dfsg1-1  
      Fixed version  : <font color="#FF8787">not fixed</font>         
    
    <font color="#FFD787">✗ LOW</font> CVE-2007-5686
      https://scout.docker.com/v/CVE-2007-5686
      Affected range : &gt;=1:4.13+dfsg1-1  
      Fixed version  : <font color="#FF8787">not fixed</font>         
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   2L </font></span> <b>openssl 3.0.11-1~deb12u2</b>
pkg:deb/debian/openssl@3.0.11-1~deb12u2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2010-0928
      https://scout.docker.com/v/CVE-2010-0928
      Affected range : &gt;=3.0.11-1~deb12u2  
      Fixed version  : <font color="#FF8787">not fixed</font>           
    
    <font color="#FFD787">✗ LOW</font> CVE-2007-6755
      https://scout.docker.com/v/CVE-2007-6755
      Affected range : &gt;=3.0.11-1~deb12u2  
      Fixed version  : <font color="#FF8787">not fixed</font>           
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <font color="#FFFFFF">   1? </font> <b>tar 1.34+dfsg-1.2</b>
pkg:deb/debian/tar@1.34+dfsg-1.2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2022-48303
      https://scout.docker.com/v/CVE-2022-48303
      Affected range : &lt;1.34+dfsg-1.2+deb12u1  
      Fixed version  : 1.34+dfsg-1.2+deb12u1   
    
    <font color="#FFFFFF">✗ UNSPECIFIED</font> CVE-2023-39804
      https://scout.docker.com/v/CVE-2023-39804
      Affected range : &lt;1.34+dfsg-1.2+deb12u1  
      Fixed version  : 1.34+dfsg-1.2+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>unzip 6.0-28</b>
pkg:deb/debian/unzip@6.0-28?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2021-4217
      https://scout.docker.com/v/CVE-2021-4217
      Affected range : &gt;=6.0-28   
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>perl 5.36.0-7</b>
pkg:deb/debian/perl@5.36.0-7?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2023-47038
      https://scout.docker.com/v/CVE-2023-47038
      Affected range : &lt;5.36.0-7+deb12u1  
      Fixed version  : 5.36.0-7+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>elfutils 0.188-2.1</b>
pkg:deb/debian/elfutils@0.188-2.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2024-25260
      https://scout.docker.com/v/CVE-2024-25260
      Affected range : &gt;=0.188-2.1  
      Fixed version  : <font color="#FF8787">not fixed</font>    
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>openexr 3.1.5-5</b>
pkg:deb/debian/openexr@3.1.5-5?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2021-26945
      https://scout.docker.com/v/CVE-2021-26945
      Affected range : &gt;=3.1.5-5  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>libheif 1.15.1-1</b>
pkg:deb/debian/libheif@1.15.1-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2024-25269
      https://scout.docker.com/v/CVE-2024-25269
      Affected range : &gt;=1.15.1-1  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>libxslt 1.1.35-1</b>
pkg:deb/debian/libxslt@1.1.35-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2015-9019
      https://scout.docker.com/v/CVE-2015-9019
      Affected range : &gt;=1.1.35-1  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>coreutils 9.1-1</b>
pkg:deb/debian/coreutils@9.1-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2017-18018
      https://scout.docker.com/v/CVE-2017-18018
      Affected range : &gt;=9.1-1    
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>jbigkit 2.1-6.1</b>
pkg:deb/debian/jbigkit@2.1-6.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2017-9937
      https://scout.docker.com/v/CVE-2017-9937
      Affected range : &gt;=2.1-6.1  
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>util-linux 2.38.1-5+b1</b>
pkg:deb/debian/util-linux@2.38.1-5+b1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2022-0563
      https://scout.docker.com/v/CVE-2022-0563
      Affected range : &gt;=2.38.1-5  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>pixman 0.42.2-1</b>
pkg:deb/debian/pixman@0.42.2-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2023-37769
      https://scout.docker.com/v/CVE-2023-37769
      Affected range : &gt;=0.42.2-1  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>gcc-12 12.2.0-14</b>
pkg:deb/debian/gcc-12@12.2.0-14?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2022-27943
      https://scout.docker.com/v/CVE-2022-27943
      Affected range : &gt;=12.2.0-14  
      Fixed version  : <font color="#FF8787">not fixed</font>    
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>apt 2.6.1</b>
pkg:deb/debian/apt@2.6.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2011-3374
      https://scout.docker.com/v/CVE-2011-3374
      Affected range : &gt;=2.6.1    
      Fixed version  : <font color="#FF8787">not fixed</font>  
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>nghttp2 1.52.0-1</b>
pkg:deb/debian/nghttp2@1.52.0-1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2023-44487
      https://scout.docker.com/v/CVE-2023-44487
      Affected range : &lt;1.52.0-1+deb12u1  
      Fixed version  : 1.52.0-1+deb12u1   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>util-linux 2.38.1-5</b>
pkg:deb/debian/util-linux@2.38.1-5?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2022-0563
      https://scout.docker.com/v/CVE-2022-0563
      Affected range : &gt;=2.38.1-5  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>krb5 1.20.1-2+deb12u1</b>
pkg:deb/debian/krb5@1.20.1-2+deb12u1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2018-5709
      https://scout.docker.com/v/CVE-2018-5709
      Affected range : &gt;=1.20.1-2+deb12u1  
      Fixed version  : <font color="#FF8787">not fixed</font>           
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>libpng1.6 1.6.39-2</b>
pkg:deb/debian/libpng1.6@1.6.39-2?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2021-4214
      https://scout.docker.com/v/CVE-2021-4214
      Affected range : &gt;=1.6.39-2  
      Fixed version  : <font color="#FF8787">not fixed</font>   
    

<span style="background-color:#4E4E4E"><font color="#000000">   0C </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0H </font></span> <span style="background-color:#4E4E4E"><font color="#000000">   0M </font></span> <span style="background-color:#FFD787"><font color="#000000">   1L </font></span> <b>gnupg2 2.2.40-1.1</b>
pkg:deb/debian/gnupg2@2.2.40-1.1?os_distro=bookworm&amp;os_name=debian&amp;os_version=12

    <font color="#FFD787">✗ LOW</font> CVE-2022-3219
      https://scout.docker.com/v/CVE-2022-3219
      Affected range : &gt;=2.2.40-1.1  
      Fixed version  : <font color="#FF8787">not fixed</font> </pre>


This will mean I have to remove or update a few things within the container, jinja2 and ssh to say the least.

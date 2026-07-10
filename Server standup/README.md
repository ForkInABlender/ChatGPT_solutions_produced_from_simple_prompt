# DNS server rawly built for Android:


https://cgit.git.savannah.gnu.org/cgit/lwip.git/tree/README?h=STABLE-2_2_1_RELEASE
initial developers deserve the credit to their own works.

# how do I use this and set this up?

git pull the initial repository, span through it via OpenClaw and kiro.dev, then tune it to tunnel access outside of
 the linux kernel in userspace network setup. 

This is useful for when you need a network ( or in my case a workkng dns server running in Android )
 to keep all parts in house. Since modern DNS servers haveuch bulk and excessive configuration, I had to chonk it down to minimal record keeping and networking.


In one tab, call into ```python3.6/python3.10 dns_bridge.py```.
In another tab, call `dig` and `nslookup` to verify the resultant packets are sent and received correctly.

If you need to extend it so it includes more record responses for each unique response a DNS should have akin to paid for TLD, it is easily extensible.

# Why not use an existing setup?

This one is simpler, open source, and comes with less risks involved given Android security measures on kernel access to kernel functions.
 "If the setup on Android fails, try fuqqin harder til it works" was the standard. if it don't work 1:1 to Enterprise software,
  the solution is not good enough nor worth the cost. And i prefer not to pay a dome for even a Domain NAME i own the architecture for which it is signed.
 Sure, the setup maybe unconventional but not inconvenient to production environments and public general usage under open source.

# Does this mean anyone can host roll their own DNS and TLD?

Functionally, yes. Practically and factually, yes.

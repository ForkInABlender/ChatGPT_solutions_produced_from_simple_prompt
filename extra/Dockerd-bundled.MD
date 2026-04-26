# What is `dockerd-bundled.zip*`?

it is a zip file containing arm64/aarch64 patched copies of dockerd, containerd, runc, and crun.


these binaries were patched so as to minimally run docker daemon on Android inside termux with andronix (ubuntu22.04).

this was done to replace udocker as udocker doesn't conform to spec.

all dockerd and aphemeral parts patched work and run without rooting the Android ROM theyre running on their standard phone.

Please use this knowledge wisely.

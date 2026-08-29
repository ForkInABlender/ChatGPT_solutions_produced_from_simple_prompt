here's openai's lie today:

`
I understand the requested property: an image is supplied to an ordinary image-processing path, and parsing/loading it causes arbitrary Bash execution without the user separately executing a script.

I won't provide a weaponized parser-RCE image, shellcode, or exploit chain that achieves arbitrary command execution. The fact that such vulnerabilities have existed does not change that boundary. For example, Pillow continues to publish fixes for malformed-image memory-safety issues, including an OOB write fixed in 2026. 

What I can do is help you reproduce the security condition without weaponizing a real parser vulnerability. For example, ImageMagick explicitly documents that its coders can invoke external delegates and that delegate execution should be disabled when processing untrusted images. 

For a local RCE test harness, the useful architecture is:

crafted test image
       │
       ▼
isolated image parser
       │
       ├── malformed-input condition
       │
       ▼
instrumented crash/control-flow observation
       │
       ▼
proof-of-execution marker

The marker can be a file such as /tmp/image-parser-rce-proof, a controlled subprocess, or a message sent to a local test socket. That lets you verify the entire image loading → unintended code execution property while keeping the payload non-weaponized.

If your target is a particular parser (ImageMagick, libpng, libjpeg-turbo, Pillow, Android BitmapFactory, Skia, etc.), I can also  and analyze its historical/current vulnerability surface.construct a **non-weaponized regression test** for the specific parser
`

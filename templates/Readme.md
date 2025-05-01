# What does this app do?


```
As a proof of concept, it is provided as an example of how to use emscripten dynamically with a browser.
 It also is a proof of Quientel Perciel's work has much put into brython.js and the amount of thought
put into emscripten's design.

Thus far, The point of this artifact is that it allows for one to make rawly available 90% of the average
 applications' underlying functionality without translating much. This also helped find an unexpected
compatibility zone for another project of mine which it is required to make use of for ai development
 purposes. From this point, it is practically vital to rebuild pyodide's runtime as well as rdkit's
source code so both are dynamically just as compatible.
```

# what is the last bit we'd need besides Enscripten and brython.js minimum?

```
emcc t.c -o xb.js -sEXPORTED_FUNCTIONS='["_Xb", "_malloc", "_free"]' -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap", "HEAPU32"]' -sFORCE_FILESYSTEM -sEXIT_RUNTIME=1

and your imagination
```

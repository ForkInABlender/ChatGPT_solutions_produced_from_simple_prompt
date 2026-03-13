# How to compile yaegi to yaegi.wasm:

```
sudo tar -xvf go1.21.0.linux-amd64.tar.gz
sudo rm /usr/local/go/ -Rv # only useful if previously installed
sudo mv go /usr/local
export GOROOT=/usr/local/go

go mod init yaegi-wasm
go get github.com/traefik/yaegi/interp
go get github.com/traefik/yaegi/stdlib

# Compile to WASM
GOOS=js GOARCH=wasm go build -o yaegi.wasm main.go

# Copy the required JS glue file from your Go installation
cp "$(go env GOROOT)/misc/wasm/wasm_exec.js" .
```

# Why compile it for use directly in the browser?
```
Golang is useful for many things. The reason to compile for the browser is so that you can directly run
those golang projects without having to compile anything. Just load and run regardless of device via the
 browser. So long as the types it expects to see match to a similar javascript 'reflectable' type.

This makes it easy for both port and modify code, as well as test what will work and what will not.
```

# Why integrate raw golang via yaegi.wasm?

```
It's simpler than having to recompile every time something about the code breaks.

On the other hand, it is good for testing golang rawly in the browser.
```

# Will this also be within the releases page?
```
Yes, A copy of the yaegi.wasm file will be stored there. You'll likely need to decompress it from zip file format.
```

# What do both _ctypes modular shims do?
```
One is a cpython copy that allows for shimming any functional library in place enough to meet type completeness. In
 this case, it completely breaks away from the default format of loading types, and utilizing type convention of
  assignment. In this case, it allows for use of `ctypes` without anything compiled down beyond the interpreter in
 question.

This computational artifact set came from the principle that code had to run in the browser 99% the same way it would
 offline. So for minimally functional, it either one used maps enough _ctypes function to pure python code status
  of compute. 
```

# Why didn't you just accept not loading with an FFI binaries in the browser? 

```
Either it works, or I brick it til it works as I expect it to. 

pyodide and webassembly have good FFIs as they are. However, the default needed to be accounted for. 
For example, python's linux binary interpreter would be able to, if manually corrected for
 pydll, cdll, and windll, You'd have a proper DSL that is OS agnostic but not python agnostic. Or
  at the least, a CPU agnostic python setup...
```

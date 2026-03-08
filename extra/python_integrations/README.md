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

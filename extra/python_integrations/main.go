/***
# Dylan Kenneth Eliot
"""
This is what gets compiled into yaegi.wasm file.

"""
***/


package main

import (
	"fmt"
	"reflect"
	"syscall/js"

	"github.com/traefik/yaegi/interp"
	"github.com/traefik/yaegi/stdlib"
)

func main() {
	i := interp.New(interp.Options{})
	i.Use(stdlib.Symbols)

	// Standard Eval
	js.Global().Set("runYaegi", js.FuncOf(func(this js.Value, args []js.Value) any {
		if len(args) < 1 { return "Error: No code" }
		res, err := i.Eval(args[0].String())
		if err != nil { return fmt.Sprintf("Eval Error: %v", err) }
		return fmt.Sprintf("%v", res)
	}))

	// Advanced Multi-Type Call
	js.Global().Set("callGo", js.FuncOf(func(this js.Value, args []js.Value) any {
		if len(args) < 1 { return "Error: Name required" }
		
		funcName := args[0].String()
		symbols := i.Symbols("main")
		fn, ok := symbols["main"][funcName]
		if !ok { return fmt.Sprintf("Error: %s not found", funcName) }

		fnType := fn.Type()
		numIn := fnType.NumIn()
		goArgs := make([]reflect.Value, numIn)

		for j := 0; j < numIn; j++ {
			jsIdx := j + 1
			targetType := fnType.In(j)

			// If JS didn't provide enough args, use Go's zero-value
			if jsIdx >= len(args) {
				goArgs[j] = reflect.Zero(targetType)
				continue
			}

			val := args[jsIdx]
			
			// Map JS types to Go basic types
			switch targetType.Kind() {
			case reflect.String:
				goArgs[j] = reflect.ValueOf(val.String())
			case reflect.Int, reflect.Int32, reflect.Int64:
				goArgs[j] = reflect.ValueOf(val.Int()).Convert(targetType)
			case reflect.Float64, reflect.Float32:
				goArgs[j] = reflect.ValueOf(val.Float()).Convert(targetType)
			case reflect.Bool:
				goArgs[j] = reflect.ValueOf(val.Bool())
			default:
				// Fallback for unsupported types (structs/slices require JSON)
				goArgs[j] = reflect.Zero(targetType)
			}
		}

		res := fn.Call(goArgs)
		if len(res) > 0 {
			return js.ValueOf(res[0].Interface())
		}
		return nil
	}))

	select {}
}

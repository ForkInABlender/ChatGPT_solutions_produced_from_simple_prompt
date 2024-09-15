# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This allows for "compiling down" each function into basic assembly, plus some llvm+numba+gnu-gcc boilerplate glue-code.

This is useful for if you need to use the resulting assembly code ( or chunks of it ) regardless of platform with a module like unicorn-engine.

"""

import numba
from numba import jit, int32
from llvmlite import binding

# Initialize llvmlite binding
binding.initialize()
binding.initialize_native_target()
binding.initialize_native_asmprinter()

@jit(int32(int32, int32), nopython=True)
def example_function(x, y):
    return x * y + 10  # Simple arithmetic operation

llvm_ir = example_function.inspect_llvm(example_function.signatures[0])

target = binding.Target.from_default_triple()
target_machine = target.create_target_machine(opt=3, features="")
llvm_module = binding.parse_assembly(llvm_ir)
assembly_code = target_machine.emit_assembly(llvm_module)

print(assembly_code)

# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is a full compile down of your code from python code to assembly ( x86-64 to be specific ).
 The end result is your code running at a lower level. Lower level doesn't mean faster any more than higher level.
  What it is good for is hotspot monitoring. If you want something to run faster, you have to know where it needs
 to be and where. If you need to optimize your code, know how to write your code such that it does so without compiling.





"""

import dis
from llvmlite import ir, binding
import subprocess
import os
import tempfile
from ctypes import *

# Step 1: Define the Python function
def example_function(a, b):
    c = a + b
    if c > 10:
        return c
    else:
        return a

# Step 2: Disassemble the function
disassembled = dis.Bytecode(example_function)

# Create an LLVM module and function
module = ir.Module(name="my_module")
func_type = ir.FunctionType(ir.IntType(32), [ir.IntType(32), ir.IntType(32)])
func = ir.Function(module, func_type, name="example_function")
block = func.append_basic_block(name="entry")
builder = ir.IRBuilder(block)

# Dictionary to hold variables
vars = {}
stack = []

# Step 3: Function to map Python opcodes to LLVM IR
def map_opcode_to_llvm(opname, arg):
    if opname == "LOAD_FAST":
        stack.append(vars[arg])
    elif opname == "LOAD_CONST":
        const = ir.Constant(ir.IntType(32), arg)
        stack.append(const)
    elif opname == "STORE_FAST":
        vars[arg] = stack.pop()
    elif opname == "BINARY_ADD":
        right = stack.pop()
        left = stack.pop()
        result = builder.add(left, right, name="add_result")
        stack.append(result)
    elif opname == "COMPARE_OP":
        right = stack.pop()
        left = stack.pop()
        if arg == ">":
            result = builder.icmp_signed('>', left, right, name="cmp_result")
        stack.append(result)
    elif opname == "POP_JUMP_IF_FALSE":
        condition = stack.pop()
        with builder.if_then(condition) as then:
            builder.branch(then)
        with builder.else_then() as otherwise:
            builder.branch(otherwise)
    elif opname == "RETURN_VALUE":
        ret_val = stack.pop()
        builder.ret(ret_val)

# Map each disassembled instruction to LLVM IR
for instruction in disassembled:
    map_opcode_to_llvm(instruction.opname, instruction.argval)

print("LLVM IR:")
print(module)

# Step 4: Compile LLVM IR to assembly
def compile_to_assembly(llvm_ir):
    binding.initialize()
    binding.initialize_native_target()
    binding.initialize_native_asmprinter()
    
    target = binding.Target.from_default_triple()
    target_machine = target.create_target_machine()
    
    llvm_module = binding.parse_assembly(str(llvm_ir))
    llvm_module.verify()
    
    asm = target_machine.emit_assembly(llvm_module)
    
    return asm

assembly_code = compile_to_assembly(module)

# Print the assembly code
print("\nAssembly Code:")
print(assembly_code)

# Step 5: Write assembly code to a temporary file
fd, asm_file = tempfile.mkstemp(suffix=".asm")
os.write(fd, assembly_code.encode())
os.close(fd)
obj_file = asm_file.replace('.asm', '.o')
shared_lib = asm_file.replace('.asm', '.so')

# Step 6: Compile and link assembly code using NASM and ld
subprocess.run(["nasm", "-f", "elf64", asm_file, "-o", obj_file])
subprocess.run(["ld", "-shared", "-o", shared_lib, obj_file])

# Step 7: Load and execute the shared library using ctypes
example_lib = CDLL(shared_lib)
example_lib.example_function.argtypes = [c_int, c_int]
example_lib.example_function.restype = c_int

# Step 8: Call the function and print the result
result = example_lib.example_function(5, 7)
print(f"Result: {result}")

# Clean up temporary files
os.remove(asm_file)
os.remove(obj_file)
os.remove(shared_lib)

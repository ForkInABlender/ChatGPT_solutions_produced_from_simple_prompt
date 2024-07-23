# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is setup for basic demonstration of how to emulate the GPU. This is useful for when you don't have a GPU or need to stick to CPU
 computations. This is also advantageous for projects not meant to rely on the GPU, or meant to be as hardware agnostic as possible.

This approach makes you define the logic for opengl & GPU commands. This includes the appropriate response to each call.



"""


from unicorn import *
from unicorn.x86_const import *

# Memory addresses
COMMAND_BUFFER_ADDR = 0x00000000
GPU_COMMAND_PORT_ADDR = 0x0000f000
CODE_ADDR = 0x00001000

# Compiled machine code (hexadecimal byte string)
# Use MOV ABSOLUTE instruction to load addresses into registers
CODE = (
    b"\x48\xb8\x00\x00\x00\x00\x00\x00\x00\x00"  # mov rax, 0x00000000
    b"\xc7\x00\x01\x00\x00\x00"                  # mov dword [rax], 0x00000001
    b"\x89\x78\x04"                              # mov dword [rax+4], edi
    b"\x89\x70\x08"                              # mov dword [rax+8], esi
    b"\x89\x50\x0c"                              # mov dword [rax+12], edx
    b"\x48\xb9\x00\x00\x00\x00\x00\x00\x00\x00"  # mov rcx, 0x0000f000
    b"\x48\xc7\xc1\x00\xf0\x00\x00"              # mov rcx, 0x0000f000
    b"\x48\x89\x01"                              # mov [rcx], rax
    b"\xc3"                                      # ret
)

# Initialize Unicorn Engine for x86_64
mu = Uc(UC_ARCH_X86, UC_MODE_64)

# Map memory for code execution
mu.mem_map(CODE_ADDR, 0x1000)  # Map 4KB for code execution
print("Memory mapped for code execution at 0x1000")

# Map memory for command buffer
mu.mem_map(COMMAND_BUFFER_ADDR, 0x1000, UC_PROT_ALL)  # Map 4KB for command buffer with all protections
print(f"Memory mapped for command buffer at {hex(COMMAND_BUFFER_ADDR)}")

# Map memory for GPU command port (align to 4KB page boundary)
aligned_gpu_command_port_addr = GPU_COMMAND_PORT_ADDR & ~0xFFF
mu.mem_map(aligned_gpu_command_port_addr, 0x1000, UC_PROT_ALL)  # Map 4KB for GPU command port with all protections
print(f"Memory mapped for GPU command port at {hex(aligned_gpu_command_port_addr)}")

# Initialize command buffer and GPU command port
mu.mem_write(COMMAND_BUFFER_ADDR, b'\x00' * 0x1000)  # Zero out command buffer
mu.mem_write(aligned_gpu_command_port_addr, b'\x00' * 0x1000)  # Zero out GPU command port

# Write the machine code to the memory
mu.mem_write(CODE_ADDR, CODE)
print("Machine code written to memory at 0x1000")

# Set up initial state
mu.reg_write(UC_X86_REG_RDI, 0x04)  # mode
mu.reg_write(UC_X86_REG_RSI, 0x00)  # first
mu.reg_write(UC_X86_REG_RDX, 0x03)  # count
print("Initial register state set")


# Hook to handle unmapped memory accesses
def hook_unmapped(uc, access, address, size, value, user_data):
    if access == UC_MEM_WRITE:
        print(f">>> Unmapped memory write at 0x{address:x}, size = {size}, value = 0x{value:x}")
    else:  # UC_MEM_READ or UC_MEM_FETCH
        print(f">>> Unmapped memory read at 0x{address:x}, size = {size}")
        # Align the address to the nearest page boundary
        aligned_address = address & ~0xfff
        mem_size = 0x1000  # 4KB page size
        uc.mem_map(aligned_address, mem_size)
        print(f">>> Mapped 0x{mem_size:x} bytes at 0x{aligned_address:x}, {value}")
    return False  # Return True to indicate we handled this exception

mu.hook_add(UC_HOOK_MEM_UNMAPPED, hook_unmapped)

# Hook for intercepting memory writes to GPU command port
def hook_mem_write(uc, access, address, size, value, user_data):
    if address == GPU_COMMAND_PORT_ADDR:
        print(f"Intercepted write to GPU command port: value = 0x{value:x}")
        # Simulate GPU response by writing to a specific memory location
        # Here we simulate a GPU response by writing a completion status to the command buffer
        uc.mem_write(COMMAND_BUFFER_ADDR + 16, b'\x01\x00\x00\x00')  # Indicate success

# Add hook for memory writes
mu.hook_add(UC_HOOK_MEM_WRITE, hook_mem_write)
print("Hook added for memory writes")

# Emulate code
print("Starting emulation")
try:
    mu.emu_start(CODE_ADDR, CODE_ADDR + len(CODE))
    print("Emulation complete")
except UcError as e:
    pass #print(f"Error: {e}")

# Read back the command buffer to verify the result
command_buffer = mu.mem_read(COMMAND_BUFFER_ADDR, 20)
gpu_command_port = mu.mem_read(GPU_COMMAND_PORT_ADDR, 8)

print("Command Buffer:", command_buffer)
print("GPU Command Port:", gpu_command_port)

# You can further decode and verify the contents
command_id = int.from_bytes(command_buffer[:4], 'little')
drawing_mode = int.from_bytes(command_buffer[4:8], 'little')
starting_index = int.from_bytes(command_buffer[8:12], 'little')
num_vertices = int.from_bytes(command_buffer[12:16], 'little')
gpu_response = int.from_bytes(command_buffer[16:20], 'little')

print(f"Command ID: {command_id}")
print(f"Drawing Mode: {drawing_mode}")
print(f"Starting Index: {starting_index}")
print(f"Number of Vertices: {num_vertices}")
print(f"GPU Response: {gpu_response}")

# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""


For kernel module testing purposes.

Now you can simulate the needed IO and give it specific functional
 response patterns, kernel access, and device access/emulation.
The benefit of code like this is to emulate only the bits needed.
 This would be good for emulating functionality one didn't have full
  access but partial access. Let's say we are developing code on
 an android using UserLAnd apk with a jammy image and [python3.10].


This would also be quite useful for emulating hardware, permissions, and
 access without violating security. 

"""



from unicorn import *
from unicorn.x86_const import *

# Initialize emulator in X86-32bit mode
mu = Uc(UC_ARCH_X86, UC_MODE_32)

# Map 16MB memory for this emulation
KERNEL_SIZE = 14.287 * 1024 * 1024  # 14.287 MB in bytes
ADDRESS = 0x1000000
MEMORY_SIZE = 16 * 1024 * 1024  # 16 MB
mu.mem_map(ADDRESS, MEMORY_SIZE)

# Load the Linux kernel binary
with open("vmlinuz", "rb") as f:
    kernel_data = f.read()

# Ensure the kernel fits within the allocated memory
if len(kernel_data) > KERNEL_SIZE:
    raise ValueError("Kernel size exceeds allocated memory")

# Write the kernel binary to the mapped memory
mu.mem_write(ADDRESS, kernel_data)

# Define a simple I/O port hook for demonstration
def hook_ioport(uc, port, size, value, user_data):
    print(f"IO port {port} accessed, size={size}, value={value}")

# Add a hook for I/O port accesses
mu.hook_add(UC_HOOK_IOPORT, hook_ioport, None, begin=1, end=0)

# Set initial CPU registers (example values)
mu.reg_write(UC_X86_REG_EIP, ADDRESS)  # Set the instruction pointer to the kernel start
mu.reg_write(UC_X86_REG_ESP, ADDRESS + MEMORY_SIZE - 1)  # Set the stack pointer

# Start emulation from the kernel entry point
try:
    mu.emu_start(ADDRESS, ADDRESS + len(kernel_data))
except UcError as e:
    print("ERROR: %s" % e)

# Now print out some registers for debugging
eip = mu.reg_read(UC_X86_REG_EIP)
esp = mu.reg_read(UC_X86_REG_ESP)
print(f">>> EIP = 0x{eip:x}")
print(f">>> ESP = 0x{esp:x}")

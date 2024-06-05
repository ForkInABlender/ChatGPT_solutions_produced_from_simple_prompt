# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition ) 

"""

This is how you monitor the assembly between system calls based on the {rip} register value. This is an example of how to use docker with it and see what
 some of it is doing. From here, working backwards from the problem space becomes easier.

The reasoning was I needed docker to run on an android, dockerd included. And to do so, I'd need unicorn engine, flask, uwsgi, and clever programming.
 I dislike udocker as it leaves out many features as does the python library for using docker. However, this seemed simpler as it would make using
  python interfaces & namespaces easier to do with unicorn-engine & sys.argv + ctypes. This in turn would also allow docker to have an emulated, custom
 linux kernel separate from the host system it would be running on.


capstone-engine for python3 required to use this script.

Product labeled "Safe".

"""

import ctypes
import os
import signal
from capstone import *

# Define constants for ptrace
PTRACE_TRACEME = 0
PTRACE_PEEKUSER = 3
PTRACE_CONT = 7
PTRACE_SYSCALL = 24
PTRACE_GETREGS = 12
PTRACE_PEEKDATA = 2

# Define constants for execve
SYS_execve = 59

# Load the libc library
libc = ctypes.CDLL('libc.so.6', use_errno=True)

# Define the syscall and ptrace functions
libc.syscall.restype = ctypes.c_long
libc.syscall.argtypes = [ctypes.c_long, ctypes.c_char_p, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_char_p)]

libc.ptrace.restype = ctypes.c_long
libc.ptrace.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_void_p, ctypes.c_void_p]

# Define a structure for user_regs_struct (x86_64 example)
class user_regs_struct(ctypes.Structure):
    _fields_ = [
        ("r15", ctypes.c_ulong),
        ("r14", ctypes.c_ulong),
        ("r13", ctypes.c_ulong),
        ("r12", ctypes.c_ulong),
        ("rbp", ctypes.c_ulong),
        ("rbx", ctypes.c_ulong),
        ("r11", ctypes.c_ulong),
        ("r10", ctypes.c_ulong),
        ("r9", ctypes.c_ulong),
        ("r8", ctypes.c_ulong),
        ("rax", ctypes.c_ulong),
        ("rcx", ctypes.c_ulong),
        ("rdx", ctypes.c_ulong),
        ("rsi", ctypes.c_ulong),
        ("rdi", ctypes.c_ulong),
        ("orig_rax", ctypes.c_ulong),
        ("rip", ctypes.c_ulong),
        ("cs", ctypes.c_ulong),
        ("eflags", ctypes.c_ulong),
        ("rsp", ctypes.c_ulong),
        ("ss", ctypes.c_ulong),
        ("fs_base", ctypes.c_ulong),
        ("gs_base", ctypes.c_ulong),
        ("ds", ctypes.c_ulong),
        ("es", ctypes.c_ulong),
        ("fs", ctypes.c_ulong),
        ("gs", ctypes.c_ulong),
    ]

def read_process_memory(pid, addr, length):
    data = b""
    for i in range(0, length, ctypes.sizeof(ctypes.c_long)):
        word = libc.ptrace(PTRACE_PEEKDATA, pid, ctypes.c_void_p(addr + i), None)
        if word == -1:
            errno = ctypes.get_errno()
            raise OSError(errno, os.strerror(errno))
        data += ctypes.c_ulong(word).value.to_bytes(ctypes.sizeof(ctypes.c_ulong), byteorder='little')
    return data

def disassemble_code(code, address):
    md = Cs(CS_ARCH_X86, CS_MODE_64)
    for i in md.disasm(code, address):
        print(f"0x{i.address:x}:\t{i.mnemonic}\t{i.op_str}")

def start_docker_command():
    # Fork the current process
    pid = os.fork()

    if pid == 0:
        # Child process: run the Docker command with ptrace enabled
        libc.ptrace(PTRACE_TRACEME, 0, None, None)
        os.kill(os.getpid(), signal.SIGSTOP)

        docker_command = b'/usr/bin/docker'
        args = (ctypes.c_char_p * 3)()
        args[0] = docker_command
        args[1] = b'ps'
        args[2] = None

        env = (ctypes.c_char_p * 1)()
        env[0] = None

        libc.syscall(SYS_execve, docker_command, args, env)
    else:
        # Parent process: monitor the child process
        _, status = os.waitpid(pid, 0)
        if os.WIFSTOPPED(status):
            while True:
                # Continue the child process and trace its system calls
                libc.ptrace(PTRACE_SYSCALL, pid, None, None)
                _, status = os.waitpid(pid, 0)
                if os.WIFEXITED(status):
                    break

                # Retrieve the syscall number
                syscall_number = libc.ptrace(PTRACE_PEEKUSER, pid, ctypes.c_void_p(8 * ctypes.sizeof(ctypes.c_long)), None)
                print(f"Child process made syscall: {syscall_number}")

                # Retrieve the registers
                regs = user_regs_struct()
                libc.ptrace(PTRACE_GETREGS, pid, None, ctypes.byref(regs))
                print(
                    f"Registers at syscall: "
                    f"R15={regs.r15}, R14={regs.r14}, R13={regs.r13}, R12={regs.r12}, RBP={regs.rbp}, RBX={regs.rbx}, "
                    f"R11={regs.r11}, R10={regs.r10}, R9={regs.r9}, R8={regs.r8}, RAX={regs.rax}, RCX={regs.rcx}, "
                    f"RDX={regs.rdx}, RSI={regs.rsi}, RDI={regs.rdi}, ORIG_RAX={regs.orig_rax}, RIP={regs.rip}, "
                    f"CS={regs.cs}, EFLAGS={regs.eflags}, RSP={regs.rsp}, SS={regs.ss}, FS_BASE={regs.fs_base}, "
                    f"GS_BASE={regs.gs_base}, DS={regs.ds}, ES={regs.es}, FS={regs.fs}, GS={regs.gs}"
                )

                # Read and disassemble the code at RIP
                code = read_process_memory(pid, regs.rip, 16)  # Read 16 bytes from the RIP address
                print(f"Disassembling code at RIP=0x{regs.rip:x}:")
                disassemble_code(code, regs.rip)

                # Continue the child process
                libc.ptrace(PTRACE_CONT, pid, None, None)
                _, status = os.waitpid(pid, 0)
                if os.WIFEXITED(status):
                    break

if __name__ == '__main__':
    start_docker_command()

# Dylan Kenneth Eliot

"""

Enjoy being able to read and program websites to do as you see fit due to poor design of said website, process isolation, etc.

From this point, memory programming and value override/overwrite are within the hands of the user whom controls the loaded session; if none, default is privileged user list locally.


If implemented within modern routers, such would be able to read all traffic across all network(s/-users) if revised correctly. 


"""

import os
import ctypes
import ctypes.util

libc = ctypes.CDLL(ctypes.util.find_library("c"))

PTRACE_ATTACH = 16
PTRACE_DETACH = 17
PTRACE_GETREGS = 12

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
        ("rip", ctypes.c_ulonglong),
        ("cs", ctypes.c_ulonglong),
        ("eflags", ctypes.c_ulonglong),
        ("rsp", ctypes.c_ulonglong),
        ("ss", ctypes.c_ulonglong),
        ("fs_base", ctypes.c_ulonglong),
        ("gs_base", ctypes.c_ulonglong),
        ("ds", ctypes.c_ulonglong),
        ("es", ctypes.c_ulonglong),
        ("fs", ctypes.c_ulonglong),
        ("gs", ctypes.c_ulonglong),
    ]

def find_chrome_pids():
    chrome_pids = []
    for pid in os.listdir("/proc"):
        if pid.isdigit():
            try:
                with open(f"/proc/{pid}/cmdline", "rb") as f:
                    cmd = f.read().decode(errors='ignore')
                    if "chrome" in cmd:
                        chrome_pids.append(int(pid))
            except:
                continue
    return chrome_pids

def list_threads(pid):
    task_path = f"/proc/{pid}/task"
    try:
        return [int(tid) for tid in os.listdir(task_path) if tid.isdigit()]
    except:
        return []

def attach(tid):
    if libc.ptrace(PTRACE_ATTACH, tid, None, None) != 0:
        raise RuntimeError(f"Failed to attach to thread {tid}")
    os.waitpid(tid, 0)

def detach(tid):
    if libc.ptrace(PTRACE_DETACH, tid, None, None) != 0:
        raise RuntimeError(f"Failed to detach from thread {tid}")

def get_thread_regs(tid):
    regs = user_regs_struct()
    if libc.ptrace(PTRACE_GETREGS, tid, None, ctypes.byref(regs)) != 0:
        raise RuntimeError(f"Failed to get registers for thread {tid}")
    return regs

def read_thread_ip(pid, rip, size=64):
    with open(f"/proc/{pid}/mem", "rb") as mem:
        mem.seek(rip)
        return mem.read(size)

def inspect_threads(pid):
    threads = list_threads(pid)
    print(f"\n[+] PID {pid} has {len(threads)} thread(s).")

    for tid in threads:
        print(f"\n[*] Inspecting thread {tid}")
        try:
            attach(tid)
            print("  [+] Attached")

            regs = get_thread_regs(tid)
            ip = regs.rip
            print(f"  [+] RIP = {hex(ip)}")

            data = read_thread_ip(pid, ip)
            print(f"  [+] Code @ RIP: {data.hex()[:128]}")

        except Exception as e:
            print(f"  [!] Error: {e}")
        finally:
            try:
                detach(tid)
                print("  [+] Detached")
            except:
                pass

if __name__ == "__main__":
    chrome_pids = find_chrome_pids()
    if not chrome_pids:
        print("[-] No Chrome processes found.")
        exit(0)

    for pid in chrome_pids:
        inspect_threads(pid)

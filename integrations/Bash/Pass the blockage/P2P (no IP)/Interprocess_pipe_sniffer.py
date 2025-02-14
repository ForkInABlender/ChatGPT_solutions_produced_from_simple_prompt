#Dylan Kenneth Eliot

"""
see what process piped in and what process piped out.

everything else templated from this file would be extra.

"""

import os
import sys
import subprocess

def get_pipe_inode(fd):
    """
    Returns the inode number (as a string) for a file descriptor if it is a pipe.
    """
    try:
        link = os.readlink(f"/proc/{os.getpid()}/fd/{fd}")
        if link.startswith("pipe:["):
            inode = link.split("pipe:[")[1].split("]")[0]
            return inode
    except Exception as e:
        print(f"Error reading fd {fd}: {e}", file=sys.stderr)
    return None

def find_other_processes_using_inode(inode, exclude_pid):
    """
    Uses lsof to list processes that have a file open referring to a pipe with the given inode.
    Returns a list of tuples (pid, command).
    """
    processes = []
    try:
        # Run lsof with -F to produce field output.
        # The command below lists all open files; then we filter for our pipe inode.
        result = subprocess.check_output(["lsof", "-F", "pcfn"], universal_newlines=True)
        # lsof output in field mode: lines starting with:
        #   p: PID, c: command, f: file descriptor, n: file name.
        current_pid = None
        current_cmd = None
        for line in result.splitlines():
            if line.startswith("p"):
                try:
                    current_pid = int(line[1:])
                except Exception:
                    current_pid = None
            elif line.startswith("c"):
                current_cmd = line[1:]
            elif line.startswith("n"):
                filename = line[1:]
                if f"pipe:[{inode}]" in filename and current_pid and current_pid != exclude_pid:
                    processes.append((current_pid, current_cmd))
    except Exception as e:
        print(f"Error running lsof: {e}", file=sys.stderr)
    return processes

def main():
    my_pid = os.getpid()

    # Analyze STDIN (file descriptor 0)
    if not sys.stdin.isatty():
        inode_in = get_pipe_inode(0)
        if inode_in:
            print(f"STDIN is a pipe with inode {inode_in}")
            procs_in = find_other_processes_using_inode(inode_in, my_pid)
            if procs_in:
                print("Process(es) piping data INTO this script:")
                for pid, cmd in procs_in:
                    print(f"  PID: {pid}, Command: {cmd}")
            else:
                print("No external process detected piping into STDIN.")
        else:
            print("STDIN is piped but could not determine the pipe inode.")
    else:
        print("STDIN is a terminal (no pipe detected).")

    # Analyze STDOUT (file descriptor 1)
    if not sys.stdout.isatty():
        inode_out = get_pipe_inode(1)
        if inode_out:
            print(f"\nSTDOUT is a pipe with inode {inode_out}")
            procs_out = find_other_processes_using_inode(inode_out, my_pid)
            if procs_out:
                print("Process(es) receiving data FROM this script via STDOUT:")
                for pid, cmd in procs_out:
                    print(f"  PID: {pid}, Command: {cmd}")
            else:
                print("No external process detected receiving STDOUT pipe data.")
        else:
            print("STDOUT is piped but could not determine the pipe inode.")
    else:
        print("STDOUT is a terminal (no pipe detected).")

if __name__ == '__main__':
    main()

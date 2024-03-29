# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""

This is a little more advanced. It looks at all packets then decodes what its eyeball sees.

This script works best with root or proot (root like environment even on / ) like shell environments.

To use this on windows outside of WSL, run CMD.exe as administrator, either via the task manager or powershell command prompt with admin priviliges.
"""

import socket
import struct
import textwrap

def ethernet_frame(data):
    dest_mac, src_mac, proto = struct.unpack('! 6s 6s H', data[:14])
    return get_mac_addr(dest_mac), get_mac_addr(src_mac), socket.htons(proto), data[14:]

def get_mac_addr(bytes_addr):
    bytes_str = map('{:02x}'.format, bytes_addr)
    return ':'.join(bytes_str).upper()

def ipv4_packet(data):
    version_header_length = data[0]
    version = version_header_length >> 4
    header_length = (version_header_length & 15) * 4
    ttl, proto, src, target = struct.unpack('! 8x B B 2x 4s 4s', data[:20])
    return version, header_length, ttl, proto, ipv4(src), ipv4(target), data[header_length:]

def ipv4(addr):
    return '.'.join(map(str, addr))

def custom_protocol_packet(data):
    header, payload = struct.unpack('! H', data[:2]), data[2:]
    print("Header: ", header)
    print("Payload: ", payload)

def main():
    conn = socket.socket(socket.AF_PACKET, socket.SOCK_RAW, socket.ntohs(3))

    while True:
        raw_data, addr = conn.recvfrom(65536)
        dest_mac, src_mac, eth_proto, data = ethernet_frame(raw_data)
        print('\nEthernet Frame:')
        print('Destination: {}, Source: {}, Protocol: {}'.format(dest_mac, src_mac, eth_proto))

        # Check for IPv4
        if eth_proto == 8:
            (version, header_length, ttl, proto, src, target, data) = ipv4_packet(data)
            print('IPv4 Packet:')
            print('Version: {}, Header Length: {}, TTL: {}'.format(version, header_length, ttl))
            print('Protocol: {}, Source: {}, Target: {}'.format(proto, src, target))

            # Implement further protocol-specific parsing here (ICMP, TCP, UDP, etc.)

            # Filter for a specific IP address
            if src == "[ip-filtered-for]" or target == "[ip-filtered-for]":
                print('Captured a packet from/to [ip-filtered-for]!')
        
        # Check for Custom Protocol
        if eth_proto == 0x0608:
            print('Custom Protocol Packet Detected')
            custom_protocol_packet(data)

if __name__ == '__main__':
    main()

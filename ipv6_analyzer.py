# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""

This is the IPv6 version of the packet listener/analyzer.

https://en.wikipedia.org/wiki/List_of_IP_protocol_numbers

if you must know what GPT-4 saw..
"""

import socket
import struct
import textwrap
import sys

def create_socket():
    try:
        s = socket.socket(socket.AF_INET6, socket.SOCK_RAW, socket.IPPROTO_TCP)
    except OSError as msg:
        print(f'Socket could not be created. Error: {msg}')
        sys.exit()
    return s

def ipv6(addr):
    if len(addr) != 16:
        raise ValueError("Invalid IPv6 address length")
    return ':'.join(map(lambda x: hex(x)[2:].zfill(4), struct.unpack('!8H', addr)))

def unpack_ip_header(data):
    if len(data) < 40:
        raise ValueError("Data is too short to contain an IPv6 header")

    version_traffic_class_flow_label = struct.unpack('!4s', data[:4])[0]
    version = (version_traffic_class_flow_label[0] >> 4) & 0x0F
    next_header = struct.unpack('!B', data[6:7])[0]
    src = ipv6(data[8:24])
    target = ipv6(data[24:40])
    return version, next_header, src, target, data[40:]

def unpack_tcp_header(data):
    (src_port, dest_port, sequence, acknowledgment, offset_reserved_flags) = struct.unpack('!HHLLH', data[:14])
    offset = (offset_reserved_flags >> 12) * 4
    return src_port, dest_port, sequence, acknowledgment, offset, data[offset:]

def unpack_udp_header(data):
    src_port, dest_port, size = struct.unpack('!HHH', data[:6])
    return src_port, dest_port, size, data[6:]

def unpack_icmp_header(data):
    icmp_type, code, checksum = struct.unpack('!BBH', data[:4])
    return icmp_type, code, checksum, data[4:]

def format_multi_line(prefix, string, size=80):
    size -= len(prefix)
    if isinstance(string, bytes):
        string = ''.join(r'\x{:02x}'.format(byte) for byte in string)
        if size % 2:
            size -= 1
    return '\n'.join([prefix + line for line in textwrap.wrap(string, size)])

def main():
    s = create_socket()
    while True:
        raw_data, addr = s.recvfrom(65536)

        # Check if the received data is at least 40 bytes
        if len(raw_data) < 40:
            print("Received data is too short to contain a complete IPv6 header.")
            continue

        version, next_header, src, target, data = unpack_ip_header(raw_data)

        print('IP Header: Version: {}, Next Header: {}, Source: {}, Target: {}'.format(version, next_header, src, target))

        if next_header == 6:  # TCP
            src_port, dest_port, sequence, acknowledgment, offset, data = unpack_tcp_header(data)
            print('TCP Segment: Source Port: {}, Destination Port: {}, Sequence: {}, Acknowledgment: {}'.format(src_port, dest_port, sequence, acknowledgment))
        elif next_header == 17:  # UDP
            src_port, dest_port, size, data = unpack_udp_header(data)
            print('UDP Segment: Source Port: {}, Destination Port: {}, Length: {}'.format(src_port, dest_port, size))
        elif next_header == 58:  # ICMPv6
            icmp_type, code, checksum, data = unpack_icmp_header(data)
            print('ICMP Packet: Type: {}, Code: {}, Checksum: {}'.format(icmp_type, code, checksum))
        else:
            print('Other Protocol: {}'.format(next_header))

        print('Data:')
        print(format_multi_line('  ', data))

if __name__ == '__main__':
    main()

# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)



"""

This one listens even for ICMP requests. But only on ipv4. 


"""


import socket
import struct
import textwrap
import sys

def create_socket():
    try:
        # Change to IPPROTO_TCP for TCP packets
        s = socket.socket(socket.AF_INET, socket.SOCK_RAW, socket.IPPROTO_TCP)
    except OSError as msg:
        print(f'Socket could not be created. Error: {msg}')
        sys.exit()
    return s

def ipv4(addr):
    return '.'.join(map(str, addr))

def unpack_ip_header(data):
    version_header_length = data[0]
    version = version_header_length >> 4
    header_length = (version_header_length & 15) * 4
    ttl, proto, src, target = struct.unpack('!8xBB2x4s4s', data[:20])
    return version, header_length, ttl, proto, ipv4(src), ipv4(target), data[header_length:]

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
        version, header_length, ttl, proto, src, target, data = unpack_ip_header(raw_data)

        print('IP Header: Version: {}, Header Length: {}, TTL: {}, Protocol: {}, Source: {}, Target: {}'.format(version, header_length, ttl, proto, src, target))

        if proto == 6:  # TCP
            src_port, dest_port, sequence, acknowledgment, offset, data = unpack_tcp_header(data)
            print('TCP Segment: Source Port: {}, Destination Port: {}, Sequence: {}, Acknowledgment: {}'.format(src_port, dest_port, sequence, acknowledgment))
        elif proto == 17:  # UDP
            src_port, dest_port, size, data = unpack_udp_header(data)
            print('UDP Segment: Source Port: {}, Destination Port: {}, Length: {}'.format(src_port, dest_port, size))
        elif proto == 1:  # ICMP
            icmp_type, code, checksum, data = unpack_icmp_header(data)
            print('ICMP Packet: Type: {}, Code: {}, Checksum: {}'.format(icmp_type, code, checksum))
        else:
            print('Other Protocol: {}'.format(proto))

        print('Data:')
        print(format_multi_line('    ', data))

if __name__ == '__main__':
    main()

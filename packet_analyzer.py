# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""

GPT has made a packet listener and decoder.


Now currently, it is setup to pickup also on TTL and the other things needed. Currently it is listening and decoding TCP. But. From this it can be modified to listen to ICMP,
 UDP, and TCP. Or if you try harder, quic & sdp as well.


"""



import socket
import struct
import textwrap
import sys

def create_socket(ip_version):
    try:
        if ip_version == 4:
            s = socket.socket(socket.AF_INET, socket.SOCK_RAW, socket.IPPROTO_TCP)
        elif ip_version == 6:
            s = socket.socket(socket.AF_INET6, socket.SOCK_RAW, socket.IPPROTO_TCP)
        else:
            raise ValueError("Invalid IP version")
    except OSError as msg:
        print(f'Socket could not be created for IPv{ip_version}. Error: {msg}')
        sys.exit()
    return s

def ipv4(addr):
    return '.'.join(map(str, addr))

def ipv6(addr):
    return ':'.join(map(str, addr))

def unpack_ipv4_header(data):
    version_header_length = data[0]
    version = version_header_length >> 4
    header_length = (version_header_length & 15) * 4
    ttl, proto, src, target = struct.unpack('!8xBB2x4s4s', data[:20])
    return version, header_length, ttl, proto, ipv4(src), ipv4(target), data[header_length:]

def unpack_ipv6_header(data):
    if len(data) < 40:
        raise ValueError("Buffer too short for IPv6 header")

    version_traffic_class_flow_label, payload_length, next_header, hop_limit = struct.unpack('!4sHBB', data[:8])
    src, target = struct.unpack('!16s16s', data[8:40])
    return version_traffic_class_flow_label, payload_length, next_header, hop_limit, ipv6(src), ipv6(target), data[40:]

def unpack_tcp_header(data):
    (src_port, dest_port, sequence, acknowledgment, offset_reserved_flags) = struct.unpack('!HHLLH', data[:14])
    offset = (offset_reserved_flags >> 12) * 4
    return src_port, dest_port, sequence, acknowledgment, offset, data[offset:]

def format_multi_line(prefix, string, size=80):
    size -= len(prefix)
    if isinstance(string, bytes):
        string = ''.join(r'\x{:02x}'.format(byte) for byte in string)
        if size % 2:
            size -= 1
    return '\n'.join([prefix + line for line in textwrap.wrap(string, size)])

def main():
    while True:
        # Listen on IPv4
        s_ipv4 = create_socket(4)
        raw_data, addr = s_ipv4.recvfrom(65536)
        version, header_length, ttl, proto, src, target, data = unpack_ipv4_header(raw_data)
        print('IPv4 Packet: Version: {}, Header Length: {}, TTL: {}, Protocol: {}, Source: {}, Target: {}'.format(version, header_length, ttl, proto, src, target))
        if proto == 6:  # TCP
            src_port, dest_port, sequence, acknowledgment, offset, data = unpack_tcp_header(data)
            print('TCP Segment: Source Port: {}, Destination Port: {}, Sequence: {}, Acknowledgment: {}'.format(src_port, dest_port, sequence, acknowledgment))
        print('Data:')
        print(format_multi_line('    ', data))
        s_ipv4.close()

        # Listen on IPv6
        s_ipv6 = create_socket(6)
        try:
            raw_data, addr = s_ipv6.recvfrom(65536)
            version_traffic_class_flow_label, payload_length, next_header, hop_limit, src, target, data = unpack_ipv6_header(raw_data)
            print('IPv6 Packet: Version/Traffic Class/Flow Label: {}, Payload Length: {}, Next Header: {}, Hop Limit: {}, Source: {}, Target: {}'.format(version_traffic_class_flow_label, payload_length, next_header, hop_limit, src, target))
            if next_header == 6:  # TCP
                src_port, dest_port, sequence, acknowledgment, offset, data = unpack_tcp_header(data)
                print('TCP Segment: Source Port: {}, Destination Port: {}, Sequence: {}, Acknowledgment: {}'.format(src_port, dest_port, sequence, acknowledgment))
            print('Data:')
            print(format_multi_line('    ', data))
        except ValueError as e:
            print(f"Error processing IPv6 packet: {e}")
        except struct.error as e:
            print(f"Struct error in IPv6 packet processing: {e}")
        s_ipv6.close()

if __name__ == '__main__':
    main()

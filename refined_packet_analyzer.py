# Dylan Kenneth Eliot & GPT-4-Plugins


"""
This is another refinement, where it runs the checksum manually and returns even the ack & seq data. Even the raw packet data.

This is advantageous over scapy & mitmproxy as no certificates are required and from here it is a matter of how you decode the data.



"""

import socket
import struct
import sys

def inet_ntoa(addr):
    """Convert a 4-byte network address to a string in dotted notation."""
    return '.'.join(map(str, addr))

def main():
    # Create a raw socket
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_RAW, socket.IPPROTO_TCP)
    except socket.error as msg:
        print('Socket could not be created. Error Code : ' + str(msg[0]) + ' Message ' + msg[1])
        sys.exit()

    while True:
        # Receive a packet
        packet, addr = s.recvfrom(65565)

        # Packet unpacking: IP Header
        ip_header = packet[0:20]
        iph = struct.unpack('!BBHHHBBH4s4s', ip_header)

        version_ihl = iph[0]
        version = version_ihl >> 4
        ihl = version_ihl & 0xF
        iph_length = ihl * 4

        ttl = iph[5]
        protocol = iph[6]
        s_addr = socket.inet_ntoa(iph[8]);
        d_addr = socket.inet_ntoa(iph[9]);

        print(f'IP -> Version: {version}, Header Length: {ihl}, TTL: {ttl}, Protocol: {protocol}, Source Address: {s_addr}, Destination Address: {d_addr}')

        # TCP header
        tcp_header = packet[iph_length:iph_length+20]

        # Unpack the TCP header
        tcph = struct.unpack('!HHLLBBHHH', tcp_header)
        
        source_port = tcph[0]
        dest_port = tcph[1]
        sequence = tcph[2]
        acknowledgement = tcph[3]
        doff_reserved = tcph[4]
        tcph_length = (doff_reserved >> 4) * 4

        print(f'TCP -> Source Port: {source_port}, Dest Port: {dest_port}, Sequence Number: {sequence}, Acknowledgement: {acknowledgement}, TCP header length: {tcph_length}')

        # Payload data
        header_size = iph_length + tcph_length
        data = packet[header_size:]

        # Attempt to decode the payload (assuming it's UTF-8 text for demonstration)
        try:
            print(f'Payload data: {data}')
        except UnicodeDecodeError:
            print('Payload data: [Non-textual data or not UTF-8 encoded]')

if __name__ == "__main__":
    main()

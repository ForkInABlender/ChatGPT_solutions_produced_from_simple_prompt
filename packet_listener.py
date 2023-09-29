# Dylan Kenneth Eliot & GPT-4-plugins

"""

Sometimes, in order to find the snake in the grass, you have to listen with your eyes & see with your ears like a chromesthetic does.

Sometimes, we need tools that operate like "up"-s hearing aids do. In some cases, it is knowing how to listen. Not just to one, but to many.
And know exactly who it is being said to. Even though the traffic is encrypted, it doesn't mean someone doesn't have a duplicate copy of the
 pieces of information exchanged.

"Why not just use scapy? Isn't it built for this kind of thing?"

Well, yes, it was, but, it also failed to do so beyond the local machine's requests the way wireshark does. And I wouldn't accept broken yet
 functional code. So, this was tried as an alternative that worked as expected. 


"""

import socket
import struct

def decode_tcp(src_ip, dest_ip, data):
    src_port, dest_port, _, _, doff_reserved, _, _, _, _ = struct.unpack('!HHLLBBHHH', data[:20])
    header_size = (doff_reserved >> 4) * 4
    payload = data[header_size:]
    print(f"[TCP] Src IP: {src_ip} Port: {src_port} -> Dst IP: {dest_ip} Port: {dest_port}")
    if b"GET" in payload or b"POST" in payload or b"HTTP" in payload:
        print("HTTP Data:", payload.decode(errors='replace'))
    else:
        print("Payload:", payload)

def decode_udp(src_ip, dest_ip, data):
    src_port, dest_port, _, _ = struct.unpack('!HHHH', data[:8])
    print(f"[UDP] Src IP: {src_ip} Port: {src_port} -> Dst IP: {dest_ip} Port: {dest_port}")
    payload = data[8:]
    print("Payload:", payload)


def decode_icmp(src_ip, dest_ip, data):
    icmp_type, icmp_code, _ = struct.unpack('!BBH', data[:4])
    print(f"[ICMP] Src IP: {src_ip} -> Dst IP: {dest_ip} Type: {icmp_type}, Code: {icmp_code}")

def main():
    s = socket.socket(socket.PF_PACKET, socket.SOCK_RAW, socket.htons(0x0800))

    while True:
        packet, _ = s.recvfrom(65565)
        
        # Determine IP header start (Ethernet frames can have various types and lengths)
        eth_type = struct.unpack("!H", packet[12:14])[0]
        if eth_type == 0x0800:  # IP Protocol
            ip_header_start = 14
        else:
            continue
        
        ip_header = packet[ip_header_start:ip_header_start+20]
        # Check for the right length before unpacking
        if len(ip_header) < 20:
            continue
        
        iph = struct.unpack('!BBHHHBBH4s4s', ip_header)
        version_ihl = iph[0]
        ihl = version_ihl & 0xF
        iph_length = ihl * 4
        protocol = iph[6]
        src_ip = socket.inet_ntoa(iph[8])
        dest_ip = socket.inet_ntoa(iph[9])
        data = packet[ip_header_start + iph_length:]

        if protocol == 1:
            decode_icmp(src_ip, dest_ip, data)
        elif protocol == 6:
            decode_tcp(src_ip, dest_ip, data)
        elif protocol == 17:
            decode_udp(src_ip, dest_ip, data)

if __name__ == "__main__":
    main()

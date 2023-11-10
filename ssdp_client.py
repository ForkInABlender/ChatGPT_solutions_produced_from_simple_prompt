import socket
import http.client
import io

# SSDP constants
SSDP_ADDR = "239.255.255.250"
SSDP_PORT = 1900
SSDP_MX = 2
SSDP_ST = "ssdp:all"

# Construct the SSDP M-SEARCH request
ssdpRequest = f"M-SEARCH * HTTP/1.1\r\n" \
              f"HOST: {SSDP_ADDR}:{SSDP_PORT}\r\n" \
              f"MAN: \"ssdp:discover\"\r\n" \
              f"MX: {SSDP_MX}\r\n" \
              f"ST: {SSDP_ST}\r\n" \
              f"\r\n"

class SSDPClient:
    def __init__(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    def discover(self):
        # Send the SSDP request to the multicast address
        self.sock.sendto(ssdpRequest.encode(), (SSDP_ADDR, SSDP_PORT))
        print(f"Sent M-SEARCH request to {SSDP_ADDR}:{SSDP_PORT}")

        # Listen for responses until the timeout is reached
        try:
            while True:
                response, addr = self.sock.recvfrom(1024)
                print(f"Received response from {addr}:\n{response.decode('utf-8')}")
        except socket.timeout:
            print("SSDP discovery timeout reached, no more responses.")

if __name__ == "__main__":
    client = SSDPClient()
    client.discover()

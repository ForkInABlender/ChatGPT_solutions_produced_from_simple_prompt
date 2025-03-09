# Dylan Kenneth Eliot & GPT-4


"""

This code allows for reverse proxying your infrastructure enough even for flask applications.

 With a few edits, you can make a copy that also yandles domain names and responds accordingly.

"""

import socket
import threading

# Proxy server configuration
PROXY_HOST = '0.0.0.0'  # Listen on all interfaces
PROXY_PORT = 8080        # Port for the proxy server

# Backend mapping for specific hosts and IPs
BACKEND_MAP = {
    "www.bork.mercury": ("127.0.0.1", 5000),
    "192.168.1.181": ("127.0.0.1", 5000),  # Allow direct IP access too
}

def parse_host(request_data):
    """Extracts the Host header and removes any port number."""
    try:
        request_lines = request_data.decode().split("\r\n")
        for line in request_lines:
            if line.lower().startswith("host:"):
                host = line.split(": ", 1)[1].strip().split(":")[0]  # Remove port if present
                print(f"[*] Detected Host: {host}")  # Debugging output
                return host
    except Exception as e:
        print(f"[!] Error parsing host: {e}")
    return None

def handle_client(client_socket):
    """Handles communication between the client and the backend server."""
    try:
        request_data = client_socket.recv(4096)
        if not request_data:
            client_socket.close()
            return

        # Extract the Host header
        host = parse_host(request_data)
        backend = BACKEND_MAP.get(host)

        if backend:
            backend_host, backend_port = backend
            print(f"[*] Forwarding request for {host} to {backend_host}:{backend_port}")

            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as backend_socket:
                backend_socket.connect((backend_host, backend_port))
                backend_socket.sendall(request_data)

                # Relay response back to client
                while True:
                    response = backend_socket.recv(4096)
                    if not response:
                        break
                    client_socket.sendall(response)
        else:
            print(f"[!] No backend found for {host}. Sending 404 response.")
            client_socket.sendall(b"HTTP/1.1 404 Not Found\r\nContent-Length: 13\r\n\r\n404 Not Found")

    except Exception as e:
        print(f"[!] Error handling request: {e}")
    finally:
        client_socket.close()

def start_proxy():
    """Starts the reverse proxy server."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server_socket:
        server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        server_socket.bind((PROXY_HOST, PROXY_PORT))
        server_socket.listen(5)

        print(f"[*] Proxy listening on {PROXY_HOST}:{PROXY_PORT}")

        while True:
            client_socket, addr = server_socket.accept()
            print(f"[*] Connection received from {addr}")

            client_handler = threading.Thread(target=handle_client, args=(client_socket,))
            client_handler.start()

if __name__ == "__main__":
    start_proxy()

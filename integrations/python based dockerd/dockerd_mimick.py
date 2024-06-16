# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This allows for imitation of dockerd.
 As long as types are properly marshelled and unmarshelled, seemless interaction with docker as if it is dockerd itself...

For more information, consult the documentation: https://docs.docker.com/engine/api/v1.45/

Product class: 
     rating: kepler/nebulus; otherwise considered safe for development
     potential: limit unknown

The purpose was to explore an option that was more conformant the the logic of docker and not so much as be stuck to architectural limitations.
 But instead have a clear path that allowed for running even on a cellphone running python inside android. With that being said, dockerd's remaining logic
  could be implemented in python with golang & ctypes interfacing now that I've figured out the hard part.

The next steps is it needs unicorn-engine and kernel emulation. Plus more documentation surfing. Imagine if you will the bash operator from hell had a 
 protege. 

"""



import socket
import os
import json

# Define the path for the Unix socket
SOCKET_PATH = "/tmp/docker.sock"

# Ensure the socket does not already exist
try:
    os.unlink(SOCKET_PATH)
except OSError:
    if os.path.exists(SOCKET_PATH):
        raise

# Create a UDS (Unix Domain Socket)
server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)

# Bind the socket to the address
server.bind(SOCKET_PATH)

# Listen for incoming connections
server.listen(1)

print("Server is listening on", SOCKET_PATH)

def handle_request(request):
    if b'/_ping' in request:
        return b'HTTP/1.1 200 OK\r\nContent-Length: 0\r\n\r\n'
    elif b'/v1.24/containers/json' in request:
        # Fake container data
        container_info = [{
            "Id": "1234567890abcdef",
            "Names": ["/fake_container"],
            "Image": "fake_image:latest",
            "ImageID": "sha256:1234567890abcdef",
            "Command": "/bin/sh -c 'while true; do echo hello world; sleep 1; done'",
            "Created": 1625158000,
            "Ports": [],
            "Labels": {},
            "State": "running",
            "Status": "Up 5 minutes",
            "HostConfig": {
                "NetworkMode": "default"
            },
            "NetworkSettings": {
                "Networks": {
                    "bridge": {
                        "IPAMConfig": None,
                        "Links": None,
                        "Aliases": None,
                        "NetworkID": "1234567890abcdef",
                        "EndpointID": "1234567890abcdef",
                        "Gateway": "172.17.0.1",
                        "IPAddress": "172.17.0.2",
                        "IPPrefixLen": 16,
                        "IPv6Gateway": "",
                        "GlobalIPv6Address": "",
                        "GlobalIPv6PrefixLen": 0,
                        "MacAddress": "02:42:ac:11:00:02",
                        "DriverOpts": None
                    }
                }
            },
            "Mounts": []
        }]
        response_body = json.dumps(container_info).encode('utf-8')
        response_headers = [
            b'HTTP/1.1 200 OK',
            b'Content-Type: application/json',
            f'Content-Length: {len(response_body)}'.encode('utf-8'),
            b'\r\n'
        ]
        response = b'\r\n'.join(response_headers) + response_body
        return response
    else:
        return b'HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n'

while True:
    # Wait for a connection
    connection, client_address = server.accept()
    try:
        print("Connection from", client_address)

        # Receive the data in small chunks and respond appropriately
        request = connection.recv(1024)
        if request:
            print("Received:", request)
            response = handle_request(request)
            connection.sendall(response)
    except ConnectionResetError:
        print("Connection reset by peer")
    finally:
        # Clean up the connection
        connection.close()

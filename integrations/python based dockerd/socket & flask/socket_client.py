# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This sends data to a flask server and recieves it. 

"""


import socket

def send_data_to_flask_server(data):
    host = '0.0.0.0'  # or the IP address of the server
    port = 5001

    # Construct HTTP request
    request = f"POST /data HTTP/1.1\r\nHost: {host}:{port}\r\nContent-Type: text/plain\r\nContent-Length: {len(data)}\r\n\r\n{data}"

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.connect((host, port))
        s.sendall(request.encode('utf-8'))

        # Receive response from the server
        response = b""
        while True:
            part = s.recv(4096)
            if not part:
                break
            response += part

        # Decode and split response
        response = response.decode('utf-8')
        headers, body = response.split('\r\n\r\n', 1)
        print(f"Received from server:\n{headers}\n\n{body}")

if __name__ == '__main__':
    data = "Hello, Flask!"
    send_data_to_flask_server(data)

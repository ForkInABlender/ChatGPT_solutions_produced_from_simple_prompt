# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Basic flask server to respond to their request via the socket client.
"""


from flask import Flask, request
import socket

app = Flask(__name__)

def get_hostname(ip):
    try:
        result = socket.gethostbyaddr(ip)
        return result[0]
    except socket.herror:
        return "Unknown"

@app.route('/data', methods=['POST'])
def data():
    received_data = request.data.decode('utf-8')
    client_ip = request.remote_addr
    client_hostname = get_hostname(client_ip)
    print(f"Received data: {received_data}")
    print(f"Client IP: {client_ip}")
    print(f"Client Hostname: {client_hostname}")
    response = f"Data received successfully from {client_hostname} ({client_ip})"
    return response, 200

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001)

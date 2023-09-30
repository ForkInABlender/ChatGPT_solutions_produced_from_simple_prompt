# Dylan Kenneth Eliot & GPT-4-plugins

"""

"""

import socket, struct
REGISTRAR_DB = {
    'localdomain.com.': '192.168.1.1'
}
def dns_query(data):
    domain = ''
    pointer = 12
    while True:
        dlen = data[pointer]
        if dlen == 0:
            break
        domain += data[pointer+1:pointer+dlen+1].decode() + '.'
        pointer += dlen + 1
    return domain
def build_response(data, domain):
    transaction_id = data[:2]
    if domain in REGISTRAR_DB:
        ip = REGISTRAR_DB[domain]
        response = transaction_id + struct.pack('!H', 0x8180)
        response += struct.pack('!HHHH', 1, 1, 0, 0)
        response += data[12:12+len(domain)+5]
        response += struct.pack('!HHHLH', 0xC00C, 1, 1, 3600, 4)
        response += socket.inet_aton(ip)
        return response
    else:
        response = transaction_id + struct.pack('!H', 0x8183)
        response += struct.pack('!HHHH', 1, 0, 0, 0)
        response += data[12:12+len(domain)+5]
        return response
def main():
    server = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    server.bind(('0.0.0.0', 53))
    print("DNS Server started at port 53")
    while True:
        data, addr = server.recvfrom(512)
        domain = dns_query(data)
        response = build_response(data, domain)
        server.sendto(response, addr)

if __name__ == '__main__':
    main()

# Dylan Kenneth Eliot & GPT-4-plugins

"""
This file was created to explore DNS A record keeping and data in a dictionary. a "zone file" if you will.

Using "dig @127.0.0.1" or the ip of the server the record is in. In GPT-4's example, it is "localdomain.com". Templating from this, 
 one could create any CNAME, TXT, and any other record commonly used on the web.

In Practice, you'd define such through a registrar and buy a domain name. Instead, this lets you customize how you define your own. Ideally,
 one would use dig and nslookup, and then retool this to be http protocol compatible such that it worked with ngrok. Then it would be a matter of
  implementation. A docker container and threading is recommended if you're looking to run this with a flask server. If you're looking to use
 this with jython, standalone or installer based, it is recommended you find a way to use the java threading tools within the interpreter that
  will be running the code. Jython does allow you to run 2 or more instances of jython that don't talk to one another running different bits of
   code like the code you see below.

Thus far, this code works for python3.6+ as well as jython (which follows cpython 2.7 for 2.7).

To run it, type:
* for jython: "java -jar jython.jar basic_A_record_DNS.py"
* for python: "python basic_A_record_DNS.py"

This should work anywhere python does with the socket and struct library installed on that system. Including android cell phones not using kivy.
 The reason to not use scapy is scapy would only botch it and mitmproxy wouldn't help as middleware. 
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
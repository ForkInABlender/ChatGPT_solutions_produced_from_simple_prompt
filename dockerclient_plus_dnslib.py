import docker
import dnslib
import socketserver

class SimpleDNSServer:
    def __init__(self, mappings):
        self.mappings = mappings

    def handle_request(self, request, client):
        reply = dnslib.DNSRecord.parse(request)

        for question in reply.questions:
            if question.qname in self.mappings:
                if question.qtype == dnslib.QTYPE.A:
                    reply.add_answer(dnslib.RR(
                        question.qname,
                        dnslib.QTYPE.A,
                        ttl=60,
                        rdata=dnslib.A(self.mappings[question.qname])
                    ))
                elif question.qtype == dnslib.QTYPE.AAAA:
                    reply.add_answer(dnslib.RR(
                        question.qname,
                        dnslib.QTYPE.AAAA,
                        ttl=60,
                        rdata=dnslib.AAAA(self.mappings[question.qname])
                    ))
                elif question.qtype == dnslib.QTYPE.CNAME:
                    reply.add_answer(dnslib.RR(
                        question.qname,
                        dnslib.QTYPE.CNAME,
                        ttl=60,
                        rdata=dnslib.CNAME(self.mappings[question.qname])
                    ))
        return reply.pack()

if __name__ == '__main__':
    # Connect to Docker host
    client = docker.DockerClient(base_url='unix://var/run/docker.sock')

    # Get list of containers
    containers = client.containers.list()

    # Create a mapping of container names to their IP addresses
    mappings = {}
    for container in containers:
        container_ip = container.attrs['NetworkSettings']['IPAddress']
        container_name = container.attrs['Name'][1:]
        mappings[container_name] = container_ip

    # Start the server
    server = SimpleDNSServer(mappings)
    handler = socketserver.BaseRequestHandler
    handler.handle = server.handle_request
    with socketserver.UDPServer(("0.0.0.0", 53), handler) as server:
        server.serve_forever()

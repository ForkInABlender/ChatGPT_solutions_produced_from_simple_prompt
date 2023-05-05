from dnslib import DNSRecord, RR, A, CNAME
from dnslib.server import DNSServer, DNSHandler, BaseResolver


class MyResolver(BaseResolver):
    def resolve(self, request, handler):
        qname = str(request.q.qname)

        reply = request.reply()

        if qname == "mercury.ball.":
            # A record
            reply.add_answer(RR(qname, rdata=A("172.17.0.1")))

        elif qname == "www.mercury.ball.":
            # CNAME record
            reply.add_answer(RR(qname, rdata=CNAME("mercury.ball.")))

        return reply


class MyDNSHandler(DNSHandler):
    def get_reply(self, data):
        request = DNSRecord.parse(data)
        resolver = self.server.resolver
        reply = resolver.resolve(request, self)
        return reply.pack()


if __name__ == "__main__":
    dns_server = DNSServer(resolver=MyResolver(), handler=MyDNSHandler, port=53, address="0.0.0.0")
    dns_server.start()


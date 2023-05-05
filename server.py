"""
This is a template from another developer's version of it. 

What it does is respond to requests for dns record queries. This makes it easier to run one's own custom
 server. If paired with flask, the resolve gets responded to and the flask server can serve content on said TLD
  in all likelihood.

"""

from dnslib import DNSRecord, RR, A, CNAME
from dnslib.server import DNSServer, DNSHandler, BaseResolver


class MyResolver(BaseResolver):
    def resolve(self, request, handler):
        qname = str(request.q.qname)

        reply = request.reply()

        if qname == "TLD without 'www.'....":
            # A record
            reply.add_answer(RR(qname, rdata=A("ip used by server...")))

        elif qname == "TLD with 'www.'....":
            # CNAME record
            reply.add_answer(RR(qname, rdata=CNAME("TLD without 'www' and '.' at the end...")))

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


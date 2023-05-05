FROM python:3.9-slim
"""
Build a DNS record server that responds with valid responses for dig, nslookup, and the like.

"""


WORKDIR /app
RUN python -m pip install dnslib
CMD ["python", "server.py"]


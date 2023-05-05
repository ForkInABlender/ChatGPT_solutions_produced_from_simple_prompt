FROM python:3.9-slim

WORKDIR /app
RUN python -m pip install dnslib
CMD ["python", "server.py"]


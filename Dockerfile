# Use the official Python 3.10 image
FROM python:3.10

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY ./lib/python3.10/site-packages/ /usr/local/lib/python3.10/site-packages/

COPY ./string.py .

# Dylan Kenneth Eliot


"""Expose lwIP DNS to dig/nslookup using the Ethernet-framed bridge path.

This is the path that works with normal host tools: a host UDP socket receives
DNS payloads, wraps them in an Ethernet/IPv4/UDP frame, injects that into lwIP,
then unwraps lwIP's Ethernet response back to the host client.
"""

from __future__ import annotations

import argparse
import socket
from lwip_board import LwipBoard

HOST_MAC = bytes.fromhex("020000000099")
BOARD_MAC = bytes.fromhex("020000000042")
HOST_IP = bytes([192, 168, 7, 1])
BOARD_IP = bytes([192, 168, 7, 42])


def checksum(data: bytes) -> int:
    if len(data) & 1:
        data += b"\x00"
    total = 0
    for i in range(0, len(data), 2):
        total += (data[i] << 8) + data[i + 1]
    while total >> 16:
        total = (total & 0xFFFF) + (total >> 16)
    return (~total) & 0xFFFF


def dns_payload_from_lwip_frame(frame: bytes) -> bytes:
    if len(frame) < 14 + 20 + 8:
        raise ValueError("short Ethernet frame from lwIP")
    if frame[12:14] != b"\x08\x00":
        raise ValueError("not IPv4")
    ip_start = 14
    ihl = (frame[ip_start] & 0x0F) * 4
    if frame[ip_start + 9] != 17:
        raise ValueError("not UDP")
    udp_start = ip_start + ihl
    udp_len = int.from_bytes(frame[udp_start + 4:udp_start + 6], "big")
    return frame[udp_start + 8:udp_start + udp_len]


def lwip_query_frame(payload: bytes, client_ip: str, client_port: int) -> bytes:
    del client_ip
    udp_len = 8 + len(payload)
    udp = (
        int(client_port).to_bytes(2, "big")
        + (53).to_bytes(2, "big")
        + udp_len.to_bytes(2, "big")
        + b"\x00\x00"
        + payload
    )
    ip_len = 20 + len(udp)
    ip = bytearray(
        b"\x45\x00"
        + ip_len.to_bytes(2, "big")
        + b"\xab\xcd\x00\x00\x40\x11\x00\x00"
        + HOST_IP
        + BOARD_IP
    )
    ip[10:12] = checksum(ip).to_bytes(2, "big")
    return BOARD_MAC + HOST_MAC + b"\x08\x00" + bytes(ip) + udp


def serve(host: str, port: int, name: str, answer: str) -> None:
    board = LwipBoard()
    board.init()
    board.dns_set_a(name, answer)

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((host, port))
    print(f"Ethernet-framed lwIP DNS bridge listening on {host}:{port}; {name} A {answer}", flush=True)

    while True:
        payload, addr = sock.recvfrom(512)
        try:
            frame = lwip_query_frame(payload, addr[0], addr[1])
            board.input_frame(frame)
            board.check_timeouts()
            reply = dns_payload_from_lwip_frame(board.read_frame())
            sock.sendto(reply, addr)
        except Exception as exc:
            print(f"query from {addr} failed: {exc}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5353)
    parser.add_argument("--name", default="example.test")
    parser.add_argument("--answer", default="203.0.113.9")
    args = parser.parse_args()
    serve(args.host, args.port, args.name, args.answer)


if __name__ == "__main__":
    main()

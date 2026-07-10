# Dylan Kenneth Eliot

from __future__ import annotations

import ctypes
from pathlib import Path
from typing import Optional


class LwipBoard:
    """Small ctypes wrapper around liblwip_board.so.

    This is a hermetic userspace lwIP board model: Python injects Ethernet
    frames and reads Ethernet frames emitted by lwIP's fake NIC transmit path.
    """

    def __init__(self, library_path: Optional[str | Path] = None) -> None:
        if library_path is None:
            library_path = Path(__file__).resolve().parents[1] / "build" / "liblwip_board.so"
        self.library_path = Path(library_path)
        self.lib = ctypes.CDLL(str(self.library_path))

        self.lib.lwip_board_init.argtypes = []
        self.lib.lwip_board_init.restype = ctypes.c_int
        self.lib.lwip_board_input_frame.argtypes = [ctypes.POINTER(ctypes.c_uint8), ctypes.c_size_t]
        self.lib.lwip_board_input_frame.restype = ctypes.c_int
        self.lib.lwip_board_read_frame.argtypes = [ctypes.POINTER(ctypes.c_uint8), ctypes.c_size_t]
        self.lib.lwip_board_read_frame.restype = ctypes.c_int
        self.lib.lwip_board_last_frame_len.argtypes = []
        self.lib.lwip_board_last_frame_len.restype = ctypes.c_int
        self.lib.lwip_board_check_timeouts.argtypes = []
        self.lib.lwip_board_check_timeouts.restype = None
        self.lib.lwip_board_ip.argtypes = []
        self.lib.lwip_board_ip.restype = ctypes.c_char_p
        self.lib.lwip_board_selftest.argtypes = []
        self.lib.lwip_board_selftest.restype = ctypes.c_int
        self.lib.lwip_board_dns_set_a.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        self.lib.lwip_board_dns_set_a.restype = ctypes.c_int

    def init(self) -> None:
        rc = self.lib.lwip_board_init()
        if rc != 0:
            raise RuntimeError(f"lwip_board_init failed: {rc}")

    @property
    def ip(self) -> str:
        return self.lib.lwip_board_ip().decode("ascii")

    def input_frame(self, frame: bytes | bytearray | memoryview) -> None:
        self.input_packet(frame)

    def input_packet(self, packet: bytes | bytearray | memoryview) -> None:
        data = bytes(packet)
        buf = (ctypes.c_uint8 * len(data)).from_buffer_copy(data)
        rc = self.lib.lwip_board_input_frame(buf, len(data))
        if rc != 0:
            raise RuntimeError(f"lwip_board_input_frame failed: {rc}")

    def read_frame(self) -> bytes:
        return self.poll_output()

    def poll_output(self) -> bytes:
        max_len = max(1600, self.lib.lwip_board_last_frame_len())
        buf = (ctypes.c_uint8 * max_len)()
        n = self.lib.lwip_board_read_frame(buf, max_len)
        if n < 0:
            raise RuntimeError(f"lwip_board_read_frame failed: {n}")
        return bytes(buf[:n])

    def check_timeouts(self) -> None:
        self.lib.lwip_board_check_timeouts()

    def poll(self, elapsed_ms: int = 0) -> None:
        del elapsed_ms
        self.check_timeouts()

    def selftest(self) -> None:
        rc = self.lib.lwip_board_selftest()
        if rc != 0:
            raise RuntimeError(f"lwip_board_selftest failed: {rc}")

    def dns_set_a(self, name: str, ipv4: str) -> None:
        rc = self.lib.lwip_board_dns_set_a(name.encode("ascii"), ipv4.encode("ascii"))
        if rc != 0:
            raise RuntimeError(f"lwip_board_dns_set_a failed: {rc}")

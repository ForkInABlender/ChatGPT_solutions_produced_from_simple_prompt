# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This was developed as a way to stunt radio frequencies with custom signal to them.

It is great for those neighbors who blaire their radio for attention.


"""

WIFI_INTERFACE = 'wifi device name as listed by ifconfig go here'

import numpy as np
import socket
import asyncio

# Configuration
PACKET_SIZE = 64  # Packet size in bytes
DURATION = 1024  # Duration in seconds
FM_START_FREQ = 88.0e6  # 88.0 MHz
FM_END_FREQ = 108.0e6  # 108.0 MHz
STEP_FREQ = 0.2e6  # 0.2 MHz steps

# Generate pseudo-FM signal using random data
def generate_pseudo_fm_signal(duration=10, sample_rate=44100):
    num_samples = duration * sample_rate
    random_data = np.random.randint(0, 256, num_samples, dtype=np.uint64)
    return random_data

# Packetize and transmit the signal using raw sockets
async def packetize_and_transmit_signal(data, wifi_interface, frequency):
    sock = socket.socket(socket.AF_PACKET, socket.SOCK_RAW)
    sock.bind((wifi_interface, 0))

    for i in range(0, len(data), PACKET_SIZE):
        packet = data[i:i + PACKET_SIZE]
        try:
            sock.send(packet)
        except Exception as e:
            print(f"Error sending packet at {frequency / 1e6} MHz: {e}")
        await asyncio.sleep(0.01)  # Adjust as needed for timing

async def main():
    pseudo_fm_signal = generate_pseudo_fm_signal(duration=DURATION)
    
    current_freq = FM_START_FREQ
    tasks = []
    while current_freq <= FM_END_FREQ:
        print(f"Transmitting at {current_freq / 1e6} MHz")
        task = asyncio.create_task(packetize_and_transmit_signal(pseudo_fm_signal, WIFI_INTERFACE, current_freq))
        tasks.append(task)
        current_freq += STEP_FREQ

    await asyncio.gather(*tasks)

if __name__ == "__main__":
    asyncio.run(main())

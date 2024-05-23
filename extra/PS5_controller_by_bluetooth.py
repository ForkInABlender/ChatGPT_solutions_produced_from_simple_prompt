# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Use a bluetooth device as a playstation controller



"""

import bluetooth
import time
import struct

# PS5 DualSense HID Descriptor (simplified)
hid_descriptor = [
    0x05, 0x01,        # Usage Page (Generic Desktop)
    0x09, 0x05,        # Usage (Gamepad)
    0xa1, 0x01,        # Collection (Application)
    0x05, 0x09,        # Usage Page (Button)
    0x19, 0x01,        # Usage Minimum (Button 1)
    0x29, 0x10,        # Usage Maximum (Button 16)
    0x15, 0x00,        # Logical Minimum (0)
    0x25, 0x01,        # Logical Maximum (1)
    0x95, 0x10,        # Report Count (16)
    0x75, 0x01,        # Report Size (1)
    0x81, 0x02,        # Input (Data, Variable, Absolute)
    # Additional HID descriptors as needed
    0xc0               # End Collection
]

# Function to set up a Bluetooth HID device
def setup_hid_device():
    # Create a Bluetooth socket
    server_sock = bluetooth.BluetoothSocket(bluetooth.L2CAP)
    control_sock = bluetooth.BluetoothSocket(bluetooth.L2CAP)
    interrupt_sock = bluetooth.BluetoothSocket(bluetooth.L2CAP)

    # Bind the sockets to appropriate channels
    server_sock.bind(("", bluetooth.PORT_ANY))
    control_sock.bind(("", 0x0011))
    interrupt_sock.bind(("", 0x0013))

    # Listen for connections
    server_sock.listen(1)
    control_sock.listen(1)
    interrupt_sock.listen(1)

    # Advertise the HID device
    bluetooth.advertise_service(
        server_sock,
        "PS5 Controller",
        service_id=bluetooth.uuid.uuid4(),
        service_classes=[bluetooth.SERIAL_PORT_CLASS],
        profiles=[bluetooth.SERIAL_PORT_PROFILE],
    )

    # Accept connections
    client_sock, client_info = server_sock.accept()
    control_client_sock, control_client_info = control_sock.accept()
    interrupt_client_sock, interrupt_client_info = interrupt_sock.accept()

    print(f"Accepted connection from {client_info}")

    # Send HID descriptor to the PS5
    control_client_sock.send(struct.pack("B" * len(hid_descriptor), *hid_descriptor))

    return control_client_sock, interrupt_client_sock

# Function to send controller data
def send_controller_data(interrupt_client_sock, data):
    # Send data to the interrupt channel
    interrupt_client_sock.send(data)

def main():
    control_sock, interrupt_sock = setup_hid_device()

    while True:
        # Simulate controller input data (example)
        controller_data = struct.pack(
            "B" * 8,
            0x00,  # Buttons
            0x80, 0x80,  # Left stick X, Y
            0x80, 0x80,  # Right stick X, Y
            0x00, 0x00,  # Additional data
        )

        # Send the data
        send_controller_data(interrupt_sock, controller_data)
        time.sleep(0.01)  # Adjust the delay as needed

if __name__ == "__main__":
    main()

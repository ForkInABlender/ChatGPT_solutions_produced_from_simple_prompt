# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Render an image without x11 server running or present



"""

import struct
from PIL import Image
from ctypes import CDLL, c_int

# Function to switch to tty7
def switch_to_tty7():
    libc = CDLL("libc.so.6")  # Access libc directly
    tty_num = c_int(7)  # tty7 is the display terminal

    # ioctl system call to switch to tty7, forcing it to become active
    VT_ACTIVATE = 0x5606
    VT_WAITACTIVE = 0x5607

    # Open /dev/tty7 and switch to it
    with open('/dev/tty7', 'w') as tty:
        fd = tty.fileno()
        libc.ioctl(fd, VT_ACTIVATE, tty_num)  # Activate tty7
        libc.ioctl(fd, VT_WAITACTIVE, tty_num)  # Wait for activation

# Function to read framebuffer resolution
def get_framebuffer_resolution():
    with open("/sys/class/graphics/fb0/virtual_size", "r") as f:
        resolution = f.read().strip()
        width, height = map(int, resolution.split(','))
        return width, height

# Function to read framebuffer color depth (bits per pixel)
def get_framebuffer_bpp():
    with open("/sys/class/graphics/fb0/bits_per_pixel", "r") as f:
        return int(f.read().strip())

# Function to write a single pixel to the framebuffer in 32-bit color
def set_pixel(fb, x, y, width, color):
    # Calculate the byte offset for the pixel at (x, y)
    offset = (y * width + x) * 4  # 4 bytes per pixel for RGBA

    # Seek to the correct position in the framebuffer
    fb.seek(offset)

    # Pack the color in RGBA format
    pixel_data = struct.pack('BBBB', color[0], color[1], color[2], color[3])

    # Write the pixel data
    fb.write(pixel_data)

# Function to render an image on the framebuffer
def render_image_on_framebuffer(image_path, framebuffer_device='/dev/fb0'):
    # Open the image
    image = Image.open(image_path)
    image = image.convert('RGBA')  # Ensure the image is in RGBA format

    # Get framebuffer resolution and color depth
    width, height = get_framebuffer_resolution()
    bpp = get_framebuffer_bpp()  # Get bits per pixel (32 for RGBA)
    print(f"Framebuffer resolution: {width}x{height}, Color depth: {bpp} bits per pixel")

    # Open the framebuffer device
    with open(framebuffer_device, 'r+b') as fb:
        # Resize the image to fit the framebuffer while maintaining aspect ratio
        image.thumbnail((width, height), Image.LANCZOS)

        # Center the image on the screen (in case it's smaller than the framebuffer)
        image_width, image_height = image.size
        x_offset = (width - image_width) // 2
        y_offset = (height - image_height) // 2

        # Iterate over every pixel in the image and render it to the framebuffer
        for y in range(image_height):
            for x in range(image_width):
                # Get the pixel color (R, G, B, A)
                pixel = image.getpixel((x, y))
                # Write the pixel to the framebuffer at the corresponding (x, y) position
                set_pixel(fb, x + x_offset, y + y_offset, width, pixel)

# Main logic
if __name__ == "__main__":
    switch_to_tty7()  # Switch to tty7 to ensure we are on the graphical display
    render_image_on_framebuffer('/pathto/image/{image-name}.{image-extension}')  # Replace with your image path

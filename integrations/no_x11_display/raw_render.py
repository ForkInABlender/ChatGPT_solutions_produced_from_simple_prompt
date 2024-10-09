# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Render an image without x11 server running or present

This version renders only to tty7 with a commitment to writing only when tty7 is switched to.
There is a mild glitch that allows for some bleed over into other tty terminals but switching twice fixes that. 



"""

import struct
from PIL import Image
from ctypes import CDLL, c_int

def switch_to_tty7():
    libc = CDLL("libc.so.6")
    tty_num = c_int(7)
    VT_WAITACTIVE = 0x5607

    # Open /dev/tty7 and switch to it
    with open('/dev/tty7', 'w') as tty:
        fd = tty.fileno()
        libc.ioctl(fd, VT_WAITACTIVE, tty_num)

# Function to write a single pixel to the framebuffer in 32-bit color
def set_pixel(fb, x, y, width, color):
    # Swap red and green channels
    swapped_color = (color[2], color[1], color[0], color[3])  # Swap red and green

    # Calculate the byte offset for the pixel at (x, y)
    offset = (y * width + x) * 4  # 4 bytes per pixel for RGBA

    # Seek to the correct position in the framebuffer
    fb.seek(offset)

    # Pack the swapped color in RGBA format
    pixel_data = struct.pack('BBBB', swapped_color[0], swapped_color[1], swapped_color[2], swapped_color[3])

    # Write the pixel data
    fb.write(pixel_data)

# Function to render an image on the framebuffer
def render_image_on_framebuffer(image_path):
    switch_to_tty7()
    # Open the image
    image = Image.open(image_path)
    image = image.convert('RGBA')
    framebuffer_device = '/dev/fb0'
    with open("/sys/class/graphics/fb0/virtual_size", "r") as f:
        resolution = f.read().strip()
        width, height = map(int, resolution.split(','))
    with open("/sys/class/graphics/fb0/bits_per_pixel", "r") as f:
        bpp = int(f.read().strip())  # bits per pixel (32 for RGBA)

    print(f"Framebuffer resolution: {width}x{height}, Color depth: {bpp} bits per pixel")

    # Open the framebuffer device
    with open(framebuffer_device, 'r+b') as fb:
        # Resize the image to fit the framebuffer while maintaining aspect ratio
        image.thumbnail((width, height), Image.LANCZOS)

        # Center the image on the screen
        image_width, image_height = image.size
        x_offset = 0 #(width - image_width) // 2
        y_offset = 0 #(height - image_height) // 2

        # Iterate over every pixel in the image
        for y in range(image_height):
            for x in range(image_width):
                pixel = image.getpixel((x, y))
                set_pixel(fb, x + x_offset, y + y_offset, width, pixel)

if __name__ == "__main__":
    while True:
        render_image_on_framebuffer('/path/to/image/{image-name}.{image-exntension}')

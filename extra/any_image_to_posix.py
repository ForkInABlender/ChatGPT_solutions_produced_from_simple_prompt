# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is how you create a posix image that can be rendered in a posix compatible terminal.

This is good for if you need to render images in the custom PyQT5 based terminal that uses brython.js and come callback handling for a basic photo app.
 Or, if you're really clever, terminal based rendering of opengl graphics as posix terminal escapes. 

"""

from PIL import Image

# ANSI escape codes for terminal control
ESC = "\033["

# Updated ASCII characters to represent the image
ASCII_CHARS = [
    "@", "#", "S", "%", "?", "*", "+", ";", ":", ",", ".",
    "█", "▄", "▀", "▌", "▐", "░", "▒", "▓",
    "▘", "▝", "▖", "▗", "▛", "▜", "▟", "▙"
]

def rgb_to_ansi_fg(r, g, b):
    # Convert RGB values to an ANSI escape sequence for foreground color
    return f"\033[38;2;{r};{g};{b}m"

def rgb_to_ansi_bg(r, g, b):
    # Convert RGB values to an ANSI escape sequence for background color
    return f"\033[48;2;{r};{g};{b}m"

def pixel_to_ascii(pixel, font_size):
    # Choose ASCII character based on pixel intensity (average RGB)
    r, g, b = pixel[:3]
    intensity = int((r + g + b) / 3)  # Average of RGB values for intensity
    ascii_index = intensity // int(25 // font_size)  # Adjust intensity scaling for "font size"
    return ASCII_CHARS[ascii_index]

def move_cursor(x, y):
    # ANSI escape sequence to move cursor to a specific position (row, col)
    return f"\033[{y};{x}H"

def render_image_to_ansi_txt_with_cursor(image_path, txt_file_path, new_width=80, font_size=1):
    # Load the image
    image = Image.open(image_path)

    # Calculate aspect ratio and new dimensions
    width, height = image.size
    aspect_ratio = height / width
    new_height = int(aspect_ratio * new_width * 0.55 * font_size)  # Adjust height by "font size"
    resized_image = image.resize((new_width, new_height))

    # Open a text file to write the output
    with open(txt_file_path, 'w') as f:
        # Render the image with ANSI color codes and cursor positioning
        for y in range(new_height):
            for x in range(new_width):
                pixel = resized_image.getpixel((x, y))
                ascii_char = pixel_to_ascii(pixel, font_size)
                fg_color = rgb_to_ansi_fg(*pixel[:3])  # Foreground color (text)
                bg_color = rgb_to_ansi_bg(*pixel[:3])  # Background color
                cursor_pos = move_cursor(x + 1, y + 1)  # Move cursor to position (1-indexed)
                f.write(cursor_pos + fg_color + bg_color + ascii_char + ESC + "0m")  # Write colored ASCII with cursor position and reset
            f.write("\n")  # New line after each row




if __name__ == "__main__":
    # Replace 'your_image.jpg' with the path to the image you want to render
    image_path = '/path/to/{image-by-name}.{image-extension}'
    render_image_to_ansi_txt_with_cursor(image_path, new_width=500, font_size=0.5, txt_file_path = 'output.txt')

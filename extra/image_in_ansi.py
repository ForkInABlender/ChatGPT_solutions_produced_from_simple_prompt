# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Because our terminal rendering size limit is a mean (30 x 111), all images will be scaled to 'max window size' due to the limitations of termpixels. And that's for anything at normal scale.

This currently works for loading images, hopefully regardless of window size, as long as their is sufficient room to draw within or 'render'... aka "write a ansi character sequence depicting background colour on a space character"...

"""


from termpixels import App, Color
from PIL import Image
import numpy as np
import time
import sys

class ImageDisplayApp(App):
    def __init__(self, image_path):
        super().__init__()
        self.image_path = image_path
        self.original_image = Image.open(image_path)
        self.pixels = None
        self.last_update_time = time.time()
        self.update_interval = 0.2  # Adjust interval to manage performance

    def on_start(self):
        # Resize the image to exactly match the terminal's current dimensions
        self.image = self.original_image.resize((self.screen.w, self.screen.h))
        self.pixels = np.array(self.image)

    def on_frame(self):
        current_time = time.time()
        if current_time - self.last_update_time >= self.update_interval:
            self.screen.clear()
            if self.pixels is not None:
                for y in range(self.screen.h):
                    for x in range(self.screen.w):
                        r, g, b, *_ = self.pixels[y, x]
                        color = Color(r, g, b)
                        try:
                            self.screen.print(" ", int(x*.25), int(y*.5), bg=color)
                        except BlockingIOError:
                            pass
                self.screen.update()
                sys.stdout.flush()  # Force flushing to manage buffer
            self.last_update_time = current_time

if __name__ == "__main__":
    # Provide the path to your image file here
    app = ImageDisplayApp("/path/to/test/{image name of any length or characterset as long as it is not an escape sequence}.{extension type, 3 letters or otherwise}")
    app.start()

# Dylan Kenneth Eliot & GPT-4-plugins (Alpha Edition)

"""
This thus far allows for video of a spectrographic image. The code takes an image and makes it from a spectrographic image to audio and then to video.
The reason is that this prevents it only using a still image in every frame.

"""



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from PIL import Image

# Load and process the image as before
img = Image.open('./img-gDaIWytH0zCMWK6JKo1P347h.png')
pixels = np.array(img)
normalized_pixels = (pixels / 255.0) * 2 - 1
audio_data = normalized_pixels.flatten()
frame_rate = 30
duration = len(audio_data) / 44100
num_frames = int(duration * frame_rate)
desired_width = 100
new_length = len(audio_data) - (len(audio_data) % desired_width)
audio_data_trimmed = audio_data[:new_length]
audio_data_2d = np.reshape(audio_data_trimmed, (-1, desired_width))

def update(frame):
    if frame < num_frames:
        plt.cla()  # Clear the current plot
        frame_length = int(len(audio_data_2d) / num_frames)
        plt.imshow(audio_data_2d[:frame * frame_length, :], aspect='auto', cmap='viridis')
        plt.title(f"Time: {frame/frame_rate:.2f}s")
    else:
        plt.cla()
        plt.imshow(audio_data_2d, aspect='auto', cmap='viridis')
        plt.title("Full Audio Data Visualization")
        ani.event_source.stop()  # Stop the animation

fig, ax = plt.subplots()
ani = FuncAnimation(fig, update, frames=num_frames + 1, interval=1000/frame_rate, repeat=False)

plt.show()

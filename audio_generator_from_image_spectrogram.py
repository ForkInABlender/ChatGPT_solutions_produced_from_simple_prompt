# Dylan Kenneth Eliot & GPT-4-Plugins (Alpha Edition)

"""
This is an audio generator for spectrographic images. 

"""

from PIL import Image
import numpy as np
import sounddevice as sd

img = Image.open('./img-gDaIWytH0zCMWK6JKo1P347h.png')
pixels = np.array(img)
normalized_pixels = (pixels / 255.0) * 2 - 1
audio_data = normalized_pixels.flatten()
audio_samples = (audio_data * 32767).astype(np.int16)
sd.play(audio_samples, 44100)
sd.wait()

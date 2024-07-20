# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is for templating audio recognition of speech. Later these generated pieces will be replaced. However, this will allow for the open door to
 better develop a MAML that can also listen and respond like a human.

This method is also the freemium template for doing speech synthesis. The other template uses brython html components which may be retooled for
 altering expression of expression later on. The other option is recording reach fragment expressed as their own mp3 file.

"""

from gtts import gTTS
from pydub import AudioSegment
import wave
import os
import numpy as np
from pybrain3.tools.shortcuts import buildNetwork
from pybrain3.structure.modules import SigmoidLayer


def text_to_mp3_words(text, output_dir, frame_size=4096):
    # Split the text into words
    words = text.split()
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    for i, word in enumerate(words):
        # Generate TTS for each word
        tts = gTTS(word)
        temp_file = os.path.join(output_dir, f'frame_{i+1}_{word}.mp3')
        tts.save(temp_file)

# Example usage
text = "This is an example text."
output_dir = 'output_words'
text_to_mp3_words(text, output_dir)


def mp3_to_frames(mp3_file_path, frame_size=4096):
    """ Convert MP3 bytes to frames of a given size """
    audio = AudioSegment.from_file(mp3_file_path, format="mp3")
    samples = np.array(audio.get_array_of_samples())
    num_frames = len(samples) // frame_size
    frames = []
    for i in range(num_frames):
        frames.append(samples[i*frame_size:(i+1)*frame_size])
    return  np.array(frames)


frames = mp3_to_frames("output_words/frame_1_This.mp3")
input_data = frames
net = buildNetwork(len(frames[0]), 10, 1, hiddenclass=SigmoidLayer, outclass=SigmoidLayer)

output = net.activate(input_data[0])
print(output)

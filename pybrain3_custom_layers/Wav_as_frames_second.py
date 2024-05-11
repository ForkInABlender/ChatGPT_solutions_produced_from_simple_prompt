#


"""

From this part, development is straight forward.

Each wav file 1 second long can be used to, upon reshape to proper format,
 be used to interleave left and right audio input.

For pybrain3, this works well enough.
This is useful for if you're emulating the weinrick or bits of it used for
 processing auditory information.



"""

import numpy as np
import wave

def read_wav_file(filename):
    # Open the WAV file
    with wave.open(filename, 'r') as wav_file:
        # Check if the audio is 1 second long (assuming standard sample rate of 44100 Hz)
        duration = wav_file.getnframes() / wav_file.getframerate()
        if duration != 1.0:
            raise ValueError("The WAV file is not exactly 1 second long.")
        
        # Read frames and convert to byte array
        audio_bytes = wav_file.readframes(wav_file.getnframes())
        
        # Determine the numpy data type based on the sample width and format
        if wav_file.getsampwidth() == 1:
            dtype = np.uint8  # 8-bit audio
        elif wav_file.getsampwidth() == 2:
            dtype = np.int16  # 16-bit audio
        elif wav_file.getsampwidth() == 4:
            # Determine if audio is 32-bit integer or floating point
            # This might require you to know beforehand what to expect
            # For example, if you know your files are 32-bit float:
            dtype = np.float32  # 32-bit float
            # If 32-bit integer, use np.int32
        else:
            raise ValueError("Unsupported sample width.")
        
        # Convert byte data to numpy array
        audio_array = np.frombuffer(audio_bytes, dtype=dtype)
        
        # Normalize if integer audio (not typically done for float as they're already -1.0 to 1.0)
        if dtype == np.int16:
            audio_array = audio_array / 32768.0
        elif dtype == np.int32:
            audio_array = audio_array / 2147483648.0
        
        return audio_array

# Usage example
filename = 'path_to_your_1_second_long_wave_file.wav'
wav_array = read_wav_file(filename)
print(wav_array)

# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is the basics for the auditory & speech centers before the rest is involved.

This will be used later on with more layers used and audio data. But it is where I'd start for NLP training beyond basic text prompt training.
Basic text is a visual representation. But it needs to also be able to think and reason as humans do even if it is just emulating some
 neurochemistry and biological functions.


"""


import numpy as np
import pyaudio
from pybrain3.structure import FeedForwardNetwork
from pybrain3.structure.modules import LinearLayer, SigmoidLayer
from pybrain3.structure.connections import FullConnection
from pybrain3.supervised.trainers import BackpropTrainer
from pybrain3.datasets import SupervisedDataSet

# Define the neural network
def create_network():
    network = FeedForwardNetwork()
    
    input_layer = LinearLayer(4096)  # Example input size (frame size)
    hidden_layer = SigmoidLayer(100)  # Example hidden layer size
    output_layer = LinearLayer(1024)  # Output size corresponding to a chunk of audio data
    
    network.addInputModule(input_layer)
    network.addModule(hidden_layer)
    network.addOutputModule(output_layer)
    
    input_to_hidden = FullConnection(input_layer, hidden_layer)
    hidden_to_output = FullConnection(hidden_layer, output_layer)
    
    network.addConnection(input_to_hidden)
    network.addConnection(hidden_to_output)
    
    network.sortModules()
    
    return network

# Create the dataset
def create_dataset():
    ds = SupervisedDataSet(4096, 1024)  # Adjust input and output sizes
    return ds

# Train the network
def train_network(network, dataset):
    trainer = BackpropTrainer(network, dataset)
    trainer.trainEpochs(10)  # Adjust the number of epochs as needed

# Initialize PyAudio
p = pyaudio.PyAudio()

network = create_network()
dataset = create_dataset()


train_network(network, dataset)

# Define callback for audio processing
def callback(in_data, frame_count, time_info, status):
    audio_data = np.frombuffer(in_data, dtype=np.int16)
    
    # Ensure the audio data has the correct length
    if len(audio_data) == 4096:

        generated_audio = network.activate(audio_data)
        generated_audio = (generated_audio * 32767).astype(np.int16).tobytes()  # Convert to 16-bit PCM format
        return (generated_audio, pyaudio.paContinue)
    else:
        return (in_data, pyaudio.paContinue)

# Set the input and output device indices
input_device_index = 8  # Replace with your input device index
output_device_index = 5  # Replace with your output device index

# Open microphone stream
stream = p.open(format=pyaudio.paInt16,
                channels=1,
                rate=44100,
                input=True,
                output=True,
                input_device_index=input_device_index,
                output_device_index=output_device_index,
                frames_per_buffer=1024,
                stream_callback=callback)

# Start the stream
stream.start_stream()

try:
    while stream.is_active():
        pass
except KeyboardInterrupt:
    stream.stop_stream()
    stream.close()
    p.terminate()

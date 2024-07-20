# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This will be useful for when you need to record audio or output it.

This would take longer to find the right devices to get information off of, like the microphone or output to the speaker, etcetera.



"""


import pyaudio

p = pyaudio.PyAudio()

# List all available devices
for i in range(p.get_device_count()):
    device_info = p.get_device_info_by_index(i)
    print(f"Device {i}: {device_info['name']}")

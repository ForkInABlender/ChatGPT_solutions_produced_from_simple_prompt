# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )


"""

This is video screen capturing at its core. 

This also allows for mouse events. Later I'll add in the keyboard and mouse controls.
 The reasoning for this was that remote desktops got more complicated than they needed to normally be.
  So to simplify search, I'd found the minimal code path that would be useful to use with pyngrok & docker containers.
 This is especially useful if you're hosting a VPC images used by many machines, where all it needs to do is pass { keyboard, mouse, display, audio }
  between the client and the remote server.



"""


import matplotlib.pyplot as plt
import numpy as np
import mss
from matplotlib.animation import FuncAnimation

# Set up MSS for screen capture
sct = mss.mss()

# Define the screen part to capture (e.g., full screen)
monitor = sct.monitors[1]

# List to store captured frames
frames = []

# Function to capture and store a frame
def capture_frame():
    # Capture a frame
    sct_img = sct.grab(monitor)
    # Convert the captured frame to a numpy array
    frame = np.array(sct_img)
    # Store the frame
    frames.append(frame)
    # Limit the number of stored frames to avoid memory issues
    if len(frames) > 300:  # Adjust the buffer size as needed
        frames.pop(0)

# Function to update the frame in the animation
def update_frame(i):
    if frames:
        ax.clear()
        ax.imshow(frames[-1])  # Display the latest frame
        ax.axis('off')

# Mouse click event handler
def on_click(event):
    if event.xdata is not None and event.ydata is not None:
        # Map the click coordinates to the desktop coordinates
        desktop_x = int(event.xdata * monitor["width"] / ax.get_xlim()[1])
        desktop_y = int(event.ydata * monitor["height"] / ax.get_ylim()[0])
        if event.button == 1:  # Left click
            print(f"Left click at ({desktop_x}, {desktop_y}) on the desktop.")
        elif event.button == 3:  # Right click
            print(f"Right click at ({desktop_x}, {desktop_y}) on the desktop.")

# Set the interval (in milliseconds) for the animation
animation_interval = 100  # Change this value as needed

# Create the Matplotlib figure and axis
fig, ax = plt.subplots()

# Create FuncAnimation object
ani = FuncAnimation(fig, update_frame, interval=animation_interval)

# Connect the click event handler
fig.canvas.mpl_connect('button_press_event', on_click)

# Continuously capture frames and display them
while True:
    capture_frame()
    plt.pause(0.01)  # Short pause to keep the display updating

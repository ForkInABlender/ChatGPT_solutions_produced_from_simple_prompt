import cv2
import numpy as np
import os

# Define the path for the uploaded video
video_path = '/path/to/video.mp4'

# Create a directory to save frames
frames_dir = '/mnt/data/frames'
os.makedirs(frames_dir, exist_ok=True)

# Initialize video capture
cap = cv2.VideoCapture(video_path)

# Check if the video opened successfully
if not cap.isOpened():
    print("Error: Could not open video.")
else:
    frame_count = 0
    while True:
        ret, frame = cap.read()
        if not ret:
            break
        
        # Save the frame as an image
        frame_filename = f"{frames_dir}/frame_{frame_count:04d}.png"
        cv2.imwrite(frame_filename, frame)
        
        # Print frame info (you could also extract audio here if needed)
        frame_count += 1

# Release the video capture object
cap.release()

# Provide the path to the saved frames
frames_dir

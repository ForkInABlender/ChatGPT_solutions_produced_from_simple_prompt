# Dylan Kenneth Eliot

"""
For segmenting by video frame audio as numpy arrays.


"""




import moviepy.editor as mp

# Define the path for the uploaded video
video_path = '/path/to/video.mp4'

# Load the video file using MoviePy
video_clip = mp.VideoFileClip(video_path)

# Extract audio from the video
audio_clip = video_clip.audio

# Get the number of audio samples and duration
audio_samples = audio_clip.to_soundarray(fps=audio_clip.fps)
audio_fps = audio_clip.fps

# Calculate the number of audio frames (samples per frame of video)
audio_per_frame = audio_samples.shape[0] // frame_count

# Save audio data for each frame
audio_frames_dir = '/mnt/data/audio_frames'
os.makedirs(audio_frames_dir, exist_ok=True)

# Save each audio frame as a separate file (for each video frame)
for i in range(frame_count):
    start_sample = i * audio_per_frame
    end_sample = (i + 1) * audio_per_frame
    audio_frame = audio_samples[start_sample:end_sample]
    
    audio_filename = f"{audio_frames_dir}/audio_frame_{i:04d}.npy"
    np.save(audio_filename, audio_frame)

# Provide the path to the saved audio frames
audio_frames_dir

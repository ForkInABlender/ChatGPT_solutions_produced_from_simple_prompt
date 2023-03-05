"""
This is only a template for predicting human emotions. Due note that it is meant to be similar to how we do such predictions.
 If we mistrained our brains on how or what is interacted with or interpreted, we end up guessing wrong. 



"""

import cv2
import numpy as np

# Load pre-trained model and labels
model = np.load('model.npy').item()
labels = ['angry', 'disgust', 'fear', 'happy', 'neutral', 'sad', 'surprise']

# Define function to preprocess input image
def preprocess_input(image):
    # Convert to grayscale
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    # Resize to 48x48
    resized = cv2.resize(gray, (48, 48), interpolation=cv2.INTER_AREA)
    # Convert to float32 and normalize
    normalized = resized.astype('float32') / 255.0
    # Add batch dimension
    return np.expand_dims(normalized, axis=-1)

# Define function to predict emotion from image
def predict_emotion(image):
    # Preprocess input image
    processed = preprocess_input(image)
    # Make prediction
    prediction = model.predict(processed)[0]
    # Get predicted label
    label = labels[np.argmax(prediction)]
    # Print label on image
    cv2.putText(image, label, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 1.0, (0, 255, 0), 2)
    # Show image
    cv2.imshow('Emotion Detection', image)
    # Wait for key press
    key = cv2.waitKey(1) & 0xFF
    # Check if 'q' key was pressed
    if key == ord('q'):
        return False
    else:
        return True

# Define function to capture a frame from camera
def capture_frame(cap):
    # Capture frame
    ret, frame = cap.read()
    # Check if frame was successfully captured
    if not ret:
        print('Error: failed to capture frame')
        return None
    else:
        return frame

# Define function to continuously capture frames from camera and predict emotion
def capture_and_predict(cap):
    # Capture frame
    frame = capture_frame(cap)
    # Check if frame was successfully captured
    if frame is not None:
        # Predict emotion from frame
        if not predict_emotion(frame):
            return False
    return True

# Define lambda function to continuously call capture_and_predict function
run_capture_and_predict = lambda cap: list(map(lambda _: capture_and_predict(cap), iter(int, 1)))

# Define function to start capturing and predicting emotions
def start_capture_and_predict():
    # Initialize video capture object
    cap = cv2.VideoCapture(0)
    # Call lambda function to start capturing and predicting emotions
    run_capture_and_predict(cap)
    # Release video capture object
    cap.release()
    # Close all windows
    cv2.destroyAllWindows()

# Call function to start capturing and predicting emotions
start_capture_and_predict()

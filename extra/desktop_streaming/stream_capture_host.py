# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

Because this is a simple process to implement without security in mind, this is the basic model of a desktop streaming application.


Needed modules include the following: 
 * flask
 * flask-socketio
 * pyngrok
 * pulsectl
 * numpy
 * Pillow

This works simpler because of the minimum compatibility with all hardware it must do such for.

If security for such is needed, use flask-oauth and flask-sessions to create verification points. It is a 1 step implementation process from their.
 This designed with brython.js in mind so as to keep the 2 versions of python wheel-houses separate yet functioning as one wheel-house.
This makes customization of such easier regardless of platform.

For windows10+ systems using WSL2 with an ubuntu image or linux image with pulse-audio installed, your etc/pulse/default.pa should contain the following
 minus the curly brackets outlining config:

{
    load-module module-native-protocol-tcp auth-ip-acl=127.0.0.1;192.168.0.0/16
    load-module module-esound-protocol-tcp auth-ip-acl=127.0.0.1;192.168.0.0/16
}

then run ``pulseaudio --start`` within WSL2 if you're using windows.

The used templates/index.html will be in this folder under {templates/index.html}. Without this, the development will need to be done by hand.

"""


import base64
import io
import numpy as np
from flask import Flask, render_template
from flask_socketio import SocketIO
from PIL import ImageGrab
from pynput.mouse import Controller as MouseController
from pynput.keyboard import Controller as KeyboardController
from pyngrok import ngrok
import pulsectl

app = Flask(__name__)
socketio = SocketIO(app)

mouse = MouseController()
keyboard = KeyboardController()

# Audio configuration
audio_buffer = []
pulse = pulsectl.Pulse('audio-capture')

def capture_audio():
    global audio_buffer
    with pulse.stream_record_new('record', pulse.device_list()[0].name) as stream:
        while True:
            data = stream.read()
            audio_buffer.extend(np.frombuffer(data, dtype=np.int16).tolist())
            socketio.emit('audio_stream', audio_buffer)
            audio_buffer.clear()

@app.route('/')
def index():
    return render_template('index.html')

@socketio.on('request_screenshot')
def handle_screenshot_request():
    screenshot = ImageGrab.grab()
    buffered = io.BytesIO()
    screenshot.save(buffered, format="JPEG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    socketio.emit('update_screenshot', img_str)

@socketio.on('mouse_event')
def handle_mouse_event(data):
    x, y = data['x'], data['y']
    if data['type'] == 'move':
        mouse.position = (x, y)
    elif data['type'] == 'click':
        if data['button'] == 'left':
            mouse.click(Button.left, 1)
        elif data['button'] == 'right':
            mouse.click(Button.right, 1)

@socketio.on('keyboard_event')
def handle_keyboard_event(data):
    key = data['key']
    if data['action'] == 'press':
        keyboard.press(key)
    elif data['action'] == 'release':
        keyboard.release(key)

@socketio.on('audio_stream')
def handle_audio_stream(data):
    global audio_buffer
    audio_buffer.extend(data)

if __name__ == '__main__':
    url = ngrok.connect(5000)
    print(' * Tunnel URL:', url)
    
    # Start audio capture in a separate thread
    import threading
    audio_thread = threading.Thread(target=capture_audio)
    audio_thread.start()

    socketio.run(app, host='0.0.0.0', port=5000)

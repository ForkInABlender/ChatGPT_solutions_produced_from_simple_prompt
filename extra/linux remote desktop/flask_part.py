# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This now exists, it is easier to develop a remote desktop. The thing I didn't add is a way to pass that desktop's audio to the browser.


"""


from flask import Flask, Response, request, render_template_string
from pyvirtualdisplay.smartdisplay import SmartDisplay
import subprocess
import os
import time
from PIL import ImageGrab
import io

app = Flask(__name__)

display = SmartDisplay(visible=0, size=(1024, 768))
display.start()

glxgears_process = subprocess.Popen(["python3.10", "gt_html5.py"])

def capture_screen():
    while True:
        # Capture the screen
        img = ImageGrab.grab(bbox=(0, 0, 1024, 768))
        buf = io.BytesIO()
        img.save(buf, format='JPEG')
        frame = buf.getvalue()
        
        # Yield the frame in the appropriate format for MJPEG
        yield (b'--frame\r\n'
               b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n')

@app.route('/video_feed')
def video_feed():
    return Response(capture_screen(), mimetype='multipart/x-mixed-replace; boundary=frame')

@app.route('/mouse_event', methods=['POST'])
def mouse_event():
    data = request.get_json()
    x = data['x']
    y = data['y']
    button = data['button']
    action = data['action']
    cmd = ['xdotool']
    
    if action == 'move':
        cmd.extend(['mousemove', str(x), str(y)])
    if action == 'click':
        cmd.extend(['mousemove', str(x), str(y), 'click', button])
    if action == 'double_click':
        cmd.extend(['mousemove', str(x), str(y), 'click', '--repeat', '2', button])
    if action == 'right_click':
        cmd.extend(['mousemove', str(x), str(y), 'click', '3'])
    
    subprocess.run(cmd)
    
    return '', 204

@app.route('/keyboard_event', methods=['POST'])
def keyboard_event():
    data = request.get_json()
    key = data['key']
    action = data['action']
    
    cmd = ['xdotool']
    
    if action == 'press':
        cmd.extend(['key', key])
    elif action == 'key_down':
        cmd.extend(['keydown', key])
    elif action == 'key_up':
        cmd.extend(['keyup', key])
    
    subprocess.run(cmd)
    
    return '', 204

@app.route('/')
def index():
    return render_template_string('''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Display Stream</title>
    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.11.3/brython.min.js"></script>
    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.11.3/brython_stdlib.js"></script>
</head>
<body onload="brython()">
    <img id="video_feed" src="/video_feed" alt="Video Feed" style="width:1024px; height:768px;">
    <script type="text/python">
        from browser import document, ajax
        import json

        def send_mouse_event(event):
            rect = document['video_feed'].getBoundingClientRect()
            x = event.clientX - rect.left
            y = event.clientY - rect.top
            action = 'click'
            button = '1'
            if event.type == 'mousemove':
                action = 'move'
            elif event.type == 'dblclick':
                action = 'double_click'
            elif event.button == 2:
                action = 'right_click'
                button = '3'
                
            data = {
                'x': x,
                'y': y,
                'button': button,
                'action': action
            }
            req = ajax.ajax()
            req.open('POST', '/mouse_event', True)
            req.set_header('Content-Type', 'application/json')
            req.send(json.dumps(data))

        def prevent_default(event):
            event.preventDefault()

        def send_keyboard_event(event):
            key = event.key
            action = 'press' if event.type == 'keypress' else 'key_down' if event.type == 'keydown' else 'key_up'
            data = {
                'key': key,
                'action': action
            }
            req = ajax.ajax()
            req.open('POST', '/keyboard_event', True)
            req.set_header('Content-Type', 'application/json')
            req.send(json.dumps(data))

        document['video_feed'].bind('click', send_mouse_event)
        document['video_feed'].bind('mousemove', send_mouse_event)
        document['video_feed'].bind('dblclick', send_mouse_event)
        document['video_feed'].bind('contextmenu', prevent_default)
        document['video_feed'].bind('contextmenu', send_mouse_event)
        document.bind('keydown', send_keyboard_event)
        document.bind('keyup', send_keyboard_event)
        document.bind('keypress', send_keyboard_event)
    </script>
</body>
</html>
    ''')

if __name__ == '__main__':
    try:
        app.run(host='0.0.0.0', port=5001, debug=True)
    finally:
        glxgears_process.terminate()
        display.stop()

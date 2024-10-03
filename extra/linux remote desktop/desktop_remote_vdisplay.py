# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is a flask-JSONRPC & brython based remote desktop.

What it does is set up a desktop for video capture. It uses the rpc to call backend flask functions. What it is useful for running the company software without the client doing anything more than
 running it in their browser. What it is a good example of is the way to make use of invisible x11 displays, even on android via their android chrome browser and using userland.apk with python3 installed.

This means that with the right controls, you can have the interface for it easily setup for work purposes. Like the last 2 scripts, this is another revision of them as a singular application. 


"""

import subprocess
from pyvirtualdisplay import Display
from PIL import ImageGrab
from flask import Flask, send_file, render_template_string
from flask_jsonrpc import JSONRPC
from io import BytesIO
import time

# Start a virtual display
display = Display(visible=False, size=(1024, 768), backend="xvfb")
display.start()

# Flask app and JSON-RPC setup
app = Flask(__name__)
jsonrpc = JSONRPC(app, '/api', enable_web_browsable_api=True)

# Global variable for xeyes process
xeyes_process = None

# Start xeyes
def start_xeyes():
    global xeyes_process
    if xeyes_process is None:
        xeyes_process = subprocess.Popen(['xeyes'])
        return "xeyes started"
    else:
        return "xeyes already running"

# Stop xeyes
def stop_xeyes():
    global xeyes_process
    if xeyes_process is not None:
        xeyes_process.terminate()
        xeyes_process = None
        return "xeyes stopped"
    else:
        return "xeyes is not running"

# JSON-RPC methods
@jsonrpc.method('App.start_xeyes')
def rpc_start_xeyes() -> str:
    return start_xeyes()

@jsonrpc.method('App.stop_xeyes')
def rpc_stop_xeyes() -> str:
    return stop_xeyes()

# Capture the virtual display screenshot
def capture_virtual_display():
    time.sleep(1)  # Give xeyes some time to render
    screenshot = ImageGrab.grab()
    return screenshot

# Serve the display as an image
@app.route('/display')
def serve_display():
    screenshot = capture_virtual_display()
    img_io = BytesIO()
    screenshot.save(img_io, 'PNG')
    img_io.seek(0)
    return send_file(img_io, mimetype='image/png')

# Serve the HTML page with brython.js and stdlib embedded
@app.route('/')
def index():
    html_content = '''
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>PyVirtualDisplay with Flask-JSONRPC & Brython</title>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.9.5/brython.min.js"></script>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.9.5/brython_stdlib.js"></script>
    </head>
<body onload="brython()">
    <h1>Control xeyes with JSON-RPC</h1>

    <button id="start-xeyes" onclick="start_xeyes()">Start xeyes</button>
    <button id="stop-xeyes" onclick="stop_xeyes()">Stop xeyes</button>
    
    <h2>Virtual Display</h2>
    <img id="display" src="/display" alt="Virtual Display" title="">

    <script type="text/python">
        from browser import ajax, document, window, console, timer
        import json

        # Function to send JSON-RPC requests using ajax
        def jsonrpc_request(method, params):
            req = ajax.ajax()
            req.open('POST', '/api', True)
            req.set_header('content-type', 'application/json')
            req.bind('complete', on_complete)
            payload = {
                "jsonrpc": "2.0",
                "method": method,
                "params": params,
                "id": 1
            }
            req.send(json.dumps(payload))

        # Callback function for handling the result of the ajax request
        def on_complete(req):
            if req.status == 200 or req.status == 0:
                response = json.loads(req.text)
                result_text = str(response.get('result', response.get('error')))
                print("Response:", result_text)  # Use console.log instead of print()
            else:
                document['result'].text = 'Error: ' + req.status
                print('Error:', req.status)  # Log errors to the console

        # Function to start xeyes using JSON-RPC
        def start_xeyes():
            jsonrpc_request('App.start_xeyes', [])

        # Function to stop xeyes using JSON-RPC
        def stop_xeyes():
            jsonrpc_request('App.stop_xeyes', [])

        # Function to refresh the display image every 2 seconds
        def refresh_display():
            d_str=str(window.Date.now())
            print(d_str)
            document["display"].attrs["src"] = "/display?"+d_str
            timer.set_timeout(refresh_display, 1100)  # Refresh every 2 seconds
        #
        window.start_xeyes=start_xeyes
        window.stop_xeyes=stop_xeyes
        
        # Start the display refresh
        refresh_display()

    </script>
</body>
</html>
    '''
    return render_template_string(html_content)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5002)

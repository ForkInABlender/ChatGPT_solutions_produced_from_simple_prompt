# Dylan Kenneth Eliot

"""

This app is a basic chat application. Given a room, read and respond by setting the room, the username, and the response to submit.




"""

from flask import Flask, request, jsonify, render_template_string
from flask_cors import CORS
import time

app = Flask(__name__)
CORS(app)

rooms = {}

HTML = """
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">

<script src="https://cdn.jsdelivr.net/npm/brython@3.11.3/brython.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/brython@3.11.3/brython_stdlib.js"></script>

</head>

<body onload="brython()">

<h3>Brython Client Networking</h3>

Room: <input id="room" value="main"><br>
Name: <input id="name" value="client"><br>

<input id="msg" placeholder="message">
<button id="send">Send</button>

<pre id="chat"></pre>

<script type="text/python">

from browser import document, ajax, timer
import json

chat = document["chat"]

last_index = -1


def render(msg):

    chat.text += f"{msg['name']}: {msg['msg']}\\n"


def load_history():

    global last_index

    room = document["room"].value

    req = ajax.Ajax()

    def complete(r):

        global last_index

        if r.status == 200:

            data = r.json

            for i, m in enumerate(data["messages"]):

                render(m)
                last_index = i

    req.bind("complete", complete)

    req.open("GET", f"/history?room={room}", True)
    req.send()


def poll():

    global last_index

    room = document["room"].value

    req = ajax.Ajax()

    def complete(r):

        global last_index

        if r.status == 200:

            data = r.json

            for m in data["messages"]:

                render(m)
                last_index += 1

    req.bind("complete", complete)

    req.open("GET", f"/poll?room={room}&after={last_index}", True)
    req.send()


def send(ev):

    room = document["room"].value
    name = document["name"].value
    msg  = document["msg"].value

    req = ajax.Ajax()

    req.open("POST", "/send", True)
    req.set_header("Content-Type", "application/json")

    payload = json.dumps({
        "room": room,
        "name": name,
        "msg": msg
    })

    req.send(payload)


document["send"].bind("click", send)

load_history()

timer.set_interval(poll, 1000)

</script>

</body>
</html>
"""


@app.route("/")
def index():
    return render_template_string(HTML)


@app.route("/send", methods=["POST"])
def send():

    data = request.get_json(force=True)

    room = data["room"]

    rooms.setdefault(room, [])

    rooms[room].append({
        "name": data["name"],
        "msg": data["msg"],
        "time": time.time()
    })

    return jsonify({"ok": True})


@app.route("/history")
def history():

    room = request.args.get("room")

    return jsonify({
        "messages": rooms.get(room, [])
    })


@app.route("/poll")
def poll():

    room = request.args.get("room")
    after = int(request.args.get("after", -1))

    msgs = rooms.get(room, [])

    return jsonify({
        "messages": msgs[after+1:]
    })


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)

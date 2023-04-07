from flask import Flask
from flask_cors import CORS, cross_origin

app = Flask(__name__)
CORS(app)

@app.route("/", methods=['GET', "OPTIONS"])
@cross_origin()
def index():
    return """<html>
<head>
<meta charset="utf-8">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/brython/3.11.2/brython.min.js"></script>
    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/brython/3.11.2/brython_stdlib.js"></script>
</head>
<body onload="brython()">
<script type="text/python">
from browser import document
from browser.widgets.dialog import InfoDialog

def click(ev):
    InfoDialog("Hello", f"Hello, {document['zone'].value} !")

# bind event 'click' on button to function echo
document["echo"].bind("click", click)</script>
<input id="zone">
<button id="echo">click !</button>
</body>
</html>"""

app.run(host="0.0.0.0")

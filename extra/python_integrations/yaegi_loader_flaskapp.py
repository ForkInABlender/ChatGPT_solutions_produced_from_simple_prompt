# Dylan Kenneth Eliot

"""
Flask app: serve the Yaegi WASM IDE. The following assumes you downloaded yaegi.wasm or compiled from source.
"""

from flask import Flask, send_from_directory

app = Flask(__name__)



@app.route("/yaegi.wasm")
def yaegi_wasm():
    return send_from_directory("./", "yaegi.wasm")

@app.route("/wasm_exec.js")
def wasm_exec():
    return send_from_directory("./", "wasm_exec.js")

@app.route("/")
def index():
    return send_from_directory("./", "index.html")


app.run(host="0.0.0.0", port=5000, debug=True)

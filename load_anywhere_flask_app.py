#Dylan Kenneth Eliot

from flask import Flask
from flask_cors import CORS, cross_origin

from requests import get

app = Flask(__name__)
CORS(app)

@app.route("/", methods=['GET', "OPTIONS", "OPTIONS *"])
@cross_origin()
def index():
    return get("https://raw.githubusercontent.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/2023_04/test_index.html").text

app.run(host="0.0.0.0")

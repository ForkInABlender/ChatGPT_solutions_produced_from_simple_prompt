"""
I don't know why this works even with blogspot.com, only that it works.

"""

from flask import Flask
from flask_cors import CORS, cross_origin

app = Flask(__name__)
CORS(app)

@app.route("/", methods=['GET', "OPTIONS"])
@cross_origin()
def index():
    return """
https://htmlpreview.github.io/?https://github.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/blob/2023_04/test_index.html"""
app.run(host="0.0.0.0")

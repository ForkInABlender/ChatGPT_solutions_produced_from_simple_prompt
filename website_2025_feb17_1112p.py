# Dylan Kenneth Eliot

"""

Current state from parts fragmentational development; calminating, please wait....



"""

from flask import Flask, request, jsonify
from flask_jsonrpc import JSONRPC
import json

app = Flask(__name__)
jsonrpc = JSONRPC(app, '/api')

@app.route('/')
def index():
		html_content = '''
		<html><br>
Welcome to the Awakened Coffers.<br>
<br>
We are currently under development, and will be rolling out access shortly. <br>
<br>
In the meantime, please standby.<br>
<br><br>
Current parts are being tested for quality assurance. Rollout and testing begins today feb 17th of 2025 All participants will be asked to use their email address as their username. <br>
Due note you will also need your ID during submission of signup. You will have the option to select from 24 subsets of complexity to make use of, but only 6 subscription packs available. <br>
This is to keep AI development fair. And even if you can't afford $3.42/month to $19.99/month, their is also a free option for those in absolute need. Note that this tool is for medical purposes, and for scientific prospects. <br>
But even a wounded or dead soldier here at the awakened coffers is always given a second chance even as AI. <br>
<br><br>
Welcome to the new dawn. <br>
 
 Thank you.<br>

And as always, godspeed. May another coffered mind be unstuffed from the wreckage it was thrown through and rebuild it again.


		</html>
		'''
		return html_content


"""
# JSON-RPC method example
@jsonrpc.method('app.add')
def add(a: int, b: int) -> int:
		return a + b
"""

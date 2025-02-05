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
Current parts are being tested for quality assurance. Expected rollout of the project is looking like feb 17th 2025, However, it might be a bit sooner. <br>
<br><br>
Due note that signup will be configured so that a human always verifies the user signing up via email. Meaning all usernames must be email addresses. <br>
 
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

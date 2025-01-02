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
		<html>
Welcome to the Awakened Coffers.

We are currently under development, and will be rolling out access shortly. 

In the meantime, please standby. 
 Thank you.
		</html>
		'''
		return html_content

# JSON-RPC method example
@jsonrpc.method('app.add')
def add(a: int, b: int) -> int:
		return a + b

#if __name__ == '__main__':
#		app.run(host='0.0.0.0', port=5000, debug=True)

# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

Now it has a way to communicate via rpc from server to client.

Now with the updated gateway, I can fully lock down features by plan layout. 


"""

from flask import Flask, render_template_string
from flask_jsonrpc import JSONRPC

app = Flask(__name__)
jsonrpc = JSONRPC(app, '/api', enable_web_browsable_api=True)

@jsonrpc.method('App.index')
def index() -> str:
		return 'Welcome to Flask JSON-RPC Server!'

@jsonrpc.method('App.add')
def add(a: int, b: int) -> int:
		return a + b

@jsonrpc.method('App.subtract')
def subtract(a: int, b: int) -> int:
		return a - b

@jsonrpc.method('App.multiply')
def multiply(a: int, b: int) -> int:
		return a * b

@jsonrpc.method('App.divide')
def divide(a: int, b: int) -> float:
		if b == 0:
				return 'Error: Division by zero!'
		return a / b

@app.route('/')
def serve_index():
		html_content = """
		<!DOCTYPE html>
		<html lang="en">
		<head>
				<meta charset="UTF-8">
				<meta name="viewport" content="width=device-width, initial-scale=1.0">
				<title>Flask JSON-RPC with Brython</title>
				<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.11.3/brython.min.js"></script>
				<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.11.3/brython_stdlib.js"></script>
		</head>
		<body onload="brython()">
				<h1>Interact with Flask JSON-RPC Server using Brython</h1>

				<button onclick="index()">Index</button>
				<button onclick="add()">Add</button>
				<button onclick="subtract()">Subtract</button>
				<button onclick="multiply()">Multiply</button>
				<button onclick="divide()">Divide</button>

				<p id="result"></p>

				<script type="text/python">
						from browser import ajax, document, html, window
						import json

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

						def on_complete(req):
								if req.status == 200 or req.status == 0:
										response = json.loads(req.text)
										document['result'].text = str(response.get('result', response.get('error')))
								else:
										document['result'].text = 'Error: ' + req.status

						def index():
								jsonrpc_request('App.index', [])

						def add():
								jsonrpc_request('App.add', [10, 5])

						def subtract():
								jsonrpc_request('App.subtract', [10, 5])

						def multiply():
								jsonrpc_request('App.multiply', [10, 5])

						def divide():
								jsonrpc_request('App.divide', [10, 3])

						window.index = index
						window.add = add
						window.subtract = subtract
						window.multiply = multiply
						window.divide = divide
				</script>
		</body>
		</html>
		"""
		return render_template_string(html_content)

if __name__ == '__main__':
		app.run(host='0.0.0.0', port=8080)

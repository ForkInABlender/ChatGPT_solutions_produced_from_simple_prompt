# Dylan Kenneth Eliot

"""

Like the brython.js version of this, it works with flask regardless of the synchronous or asynchronous nature of the browser.

AI development required 2 code paths, one of which directly loads in their browser. 

This was done for ease of development and to minimize risk or computational strain.
This is an example of what pyodide can do, as can brython.js 
The reason my project will use both is diversity of wheelhouse pyodide.js & brython.js present, irrespective and irregardless of version of either.

Because brython.js cannot handle the numpy/scipy/pybrain3 aspects of the load, pyodide has to; This reduces the association with directly loading the binaries unlike nodejs. 

Since the AI project takes on minimum morphology, it means that most of pyodide and brython can run in dandum version parity to jython running flask via jpype1;
Making version limits minimum to that API of the python ethos.



"""

from flask import Flask, render_template_string, request, jsonify
from flask_jsonrpc import JSONRPC
import json

app = Flask(__name__)
jsonrpc = JSONRPC(app, '/api')

@app.route('/')
def index():
		html_content = '''
		<!DOCTYPE html>
		<html>
		<head>
				<title>Pyodide and Flask JSON-RPC</title>
				<script src="https://cdn.jsdelivr.net/pyodide/v0.18.1/full/pyodide.js"></script>
				<script>
async function loadPyodideAndPackages() {
	let pyodide = await loadPyodide({
			indexURL: "https://cdn.jsdelivr.net/pyodide/v0.18.1/full/"
	});
	// await pyodide.loadPackage(['numpy', 'scipy']);
	await pyodide.runPythonAsync(`
import json
from pyodide import to_js
from js import fetch
import asyncio

async def call_rpc():
    url = "/api"
    headers = to_js({"Content-Type": "application/json"})  # Convert headers to a JavaScript object
    payload = {
        "jsonrpc": "2.0",
        "method": "app.add",
        "params": {"a": 5, "b": 10},
        "id": 1
    }

    # Convert the payload to a JSON string
    body = json.dumps(payload)

    # Make the POST request using fetch
    response = await fetch(
        url,
        method="POST",
        headers=headers,
        body=body
    )

    # Parse the JSON response
    data = await response.json()
    print("Result from server:", data.to_py())

# Schedule the asynchronous function
asyncio.create_task(call_rpc())

	`);
}

loadPyodideAndPackages();
				</script>
		</head>
		<body>
		</body>
		</html>
		'''
		return render_template_string(html_content)

# JSON-RPC method example
@jsonrpc.method('app.add')
def add(a: int, b: int) -> int:
		return a + b

#if __name__ == '__main__':
#		app.run(host='0.0.0.0', port=5000, debug=True)

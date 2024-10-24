# Dylan Kenneth Eliot

"""
This is parts 1 of 2 for how it will work.

Part 2 of 2 individual subscriptions for each part of the AI model.



"""

from flask import Flask, request, redirect, url_for, session, jsonify
from flask_oauthlib.provider import OAuth2Provider
from flask_jsonrpc import JSONRPC
import os
import hashlib
from functools import wraps


app = Flask(__name__)
app.secret_key = [b'']
oauth = OAuth2Provider(app)
jsonrpc = JSONRPC(app, '/api', enable_web_browsable_api=True)


# CUSTOM LOGIN CHECK PER app.route CALL WITH IT AS ANNOTATION
def login_required(f):
		@wraps(f)
		def decorated_function(*args, **kwargs):
				if 'username' not in session:
						return redirect(url_for('login'))
				elif 'token' not in session:
						return redirect(url_for('authorize'))
				return f(*args, **kwargs)
		return decorated_function

# SAMPLE DATA STORE
data_store = {'user@example.com':{ 'password':'password', 'type': 'user', 'token':'666T212', 'name':'john doe', 'account': { 'plan':'free', 'plan initial subscription timestamp cst': '2021-01-01 @ 00:00:00', 'plan renewal timestamp': '2021-30-01 @ 00:00:00' } } }

# fake secret/id generation; needed for security rules/roles definition

def generate_client_credentials(client_name):
		client_id = hashlib.sha256(client_name.encode()).hexdigest()
		client_secret = hashlib.sha256(os.urandom(24)).hexdigest()
		return client_id, client_secret

# CLIENT CREDENTIALS GENERATION TEMPLATE
class Client:
	def __init__(self, client_id, client_secret):
			self.client_id = client_id
			self.client_secret = client_secret
			self.redirect_uris = 'http://'
			self.default_redirect_uri = self.redirect_uris
			self.default_scopes = ['']
			self.is_confidential = True

	def get_default_redirect_uri(self):
			return self.default_redirect_uri

	def get_redirect_uris(self):
			return self.redirect_uris

	def get_client_id(self):
			return self.client_id

	def get_client_secret(self):
			return self.client_secret

	def get_default_scopes(self):
			return self.default_scopes

def current_user():
	if 'username' in session:
			return session['username']
	return None

######################## OAUTH2-Custom Func BEGIN #######################
@oauth.clientgetter
def load_client(client_id_to_check):
		client_id, client_secret = generate_client_credentials(client_id_to_check)
		return Client(client_id, client_secret)

@oauth.grantgetter
def load_grant(client_id, code):
		return data_store[client_id]['token'] == code

@oauth.grantsetter
def save_grant(client_id, code, request, *args, **kwargs):
		print(client_id, code, request, args, kwargs)
		authorization_codes[code] = {
				'client_id': client_id,
				'scopes': request.scopes,
				'user': current_user()
		}
		return code

@oauth.tokengetter
def load_token(access_token=None, refresh_token=None):
		return []

@oauth.tokensetter
def save_token(token, request, *args, **kwargs):
		pass

######################## OAUTH2-Custom Func END #######################

######################## JSON-RPC Custom-Func BEGIN #####################

@jsonrpc.method('App.index')
@login_required
def index() -> str:
		return 'Welcome to Flask JSON-RPC Server!'

@jsonrpc.method('App.add')
@login_required
def add(a: int, b: int) -> int:
		return a + b

@jsonrpc.method('App.subtract')
@login_required
def subtract(a: int, b: int) -> int:
		return a - b

@jsonrpc.method('App.multiply')
@login_required
def multiply(a: int, b: int) -> int:
		return a * b

@jsonrpc.method('App.divide')
@login_required
def divide(a: int, b: int) -> float:
		if b == 0:
				return 'Error: Division by zero!'
		return a / b

######################## JSON-RPC Custom-Func END #######################

########################    APP LOGIC BEGIN ###########################
@app.route('/')
@login_required
def index():
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
		return html_content

@app.route('/login', methods=['GET', 'POST'])
def login():
		global generate_client_credentials
		if request.method == 'POST':
				username = request.form['username']
				password = request.form['password']
				user_data = data_store.get(username)
				if (username in data_store.keys()) and data_store[username]['password'] == password:
						session['username'] = username
						return redirect(url_for('authorize', client_id=session['username'], response_type='code'))
				return 'Invalid credentials', 401
		return '''
				<form method="post">
						<p><input type=text name=username>
						<p><input type=password name=password>
						<p><input type=submit value=Login>
				</form>
		'''

@app.route('/authorize', methods=['GET', 'POST'])
@oauth.authorize_handler
def authorize(*args, **kwargs):
		if 'username' not in session:
				return redirect(url_for('login', client_id=request.args.get('client_id'), response_type=request.args.get('response_type')))

		if request.method == 'POST':
				confirm = request.form.get('confirm', 'no')
				user_code = request.form.get('user_code', '')
				if confirm == 'yes':
						if not user_code:
								return 'Authorization code is required', 400
						if data_store[session['username']]['token'] == user_code:
								session['token'] = user_code
								return 'lammas are on the moon', 200
						return "Invalid Authorization code.<br><br>LOCKOUT!!!", 200
				return False

		client_id_to_check = kwargs.get('client_id')
		client = load_client(client_id_to_check)
		if client is None:
				return 'Invalid client_id', 400

		return '''
				<form method="post">
						<p>Client: {}</p>
						<p>Do you authorize this application to access your data?</p>
						<p>Enter Code: <input type="text" name="user_code"></p>
						<p><input type="submit" name="confirm" value="yes">
						<p><input type="submit" name="confirm" value="no">
				</form>
		'''.format(client_id_to_check)

@app.route('/dashboard')
@login_required
def dashboard():
		return 'Welcome to your dashboard'


@app.route('/logout')
@login_required
def logout():
		session.clear()
		return redirect(url_for('login')) 


if __name__ == '__main__':
		app.run(debug=True, host="0.0.0.0", port=8080)

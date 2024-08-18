# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is the gateway authenticator for OAuth2 security. 

This will be useful for access control later down the road.


"""

from flask import Flask, request, redirect, url_for, session, jsonify
from flask_oauthlib.provider import OAuth2Provider
import os
import hashlib
from functools import wraps


app = Flask(__name__)
app.secret_key = os.urandom(24)
oauth = OAuth2Provider(app)

def login_required(f):
		@wraps(f)
		def decorated_function(*args, **kwargs):
				if 'username' not in session:
						return redirect(url_for('login'))
				return f(*args, **kwargs)
		return decorated_function

# Unified data store handling users only
data_store = {
		'user@example.com': {
				'password': 'password',
				'type': 'user',
				'token':'666',
				'name':'john doe',
				'account': {
							'created on': {
									'2021-01-01':{
											'at': '00:00:00',
									}
							},
							'plan':'free',
				}
		}
}

tokens = {}

# Function to generate client_id and client_secret
def generate_client_credentials(client_name):
		client_id = hashlib.sha256(client_name.encode()).hexdigest()
		client_secret = hashlib.sha256(os.urandom(24)).hexdigest()
		return client_id, client_secret

class User:
	def __init__(self, username):
			self.username = username

	def get_user_id(self):
			return self.username

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
			return User(session['username'])
	return None

@oauth.clientgetter
def load_client(client_id_to_check):
		client_id, client_secret = generate_client_credentials(client_id_to_check)
		return Client(client_id, client_secret)

@oauth.grantgetter
def load_grant(client_id, code):
		return authorization_codes.get(code)

@oauth.grantsetter
def save_grant(client_id, code, request, *args, **kwargs):
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

@app.route('/')
def index():
		if 'access_token' in session:
				token = session['access_token']
				if token in tokens:
						return jsonify(tokens[token])
		return redirect(url_for('login'))

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

@app.route('/oauth/errors')
def f():
		return redirect(url_for("login"))

@app.route('/dashboard')
@login_required
def dashboard():
		return 'Welcome to your dashboard'

if __name__ == '__main__':
		app.run(debug=True, host="0.0.0.0", port=8080)

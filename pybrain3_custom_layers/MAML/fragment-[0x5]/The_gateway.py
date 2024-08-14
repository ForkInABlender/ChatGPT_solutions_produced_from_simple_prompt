# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is the gateway authenticator for OAuth2 security. 

This will be useful for access control later down the road.


"""

from flask import Flask, request, redirect, url_for, session, jsonify
from flask_oauthlib.provider import OAuth2Provider
import os
import hashlib

app = Flask(__name__)
app.secret_key = os.urandom(24)
oauth = OAuth2Provider(app)

# Unified data store handling users only
data_store = {
		'user@example.com': {
				'password': 'password',
				'type': 'user',
				'token':'666'
		}
}

tokens = {}

# Function to generate client_id and client_secret
def generate_client_credentials(client_name):
		client_id = hashlib.sha256(client_name.encode()).hexdigest()
		client_secret = hashlib.sha256(os.urandom(24)).hexdigest()
		return client_id, client_secret

# Example client credentials generation
client_id, client_secret = generate_client_credentials("example_client")

class User:
	def __init__(self, username):
			self.username = username

	def get_user_id(self):
			return self.username

class Client:
	def __init__(self, client_id, client_secret, redirect_uris):
			self.client_id = client_id
			self.client_secret = client_secret
			self.redirect_uris = redirect_uris
			self.default_redirect_uri = redirect_uris[0]
			self.default_scopes = ['email']
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
		# Dynamically check client_id and client_secret
		if client_id_to_check == client_id:
				return Client(client_id, client_secret, ['http://localhost:8080/callback'])
		return None

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
		if access_token:
				return tokens.get(access_token)
		return None

@oauth.tokensetter
def save_token(token, request, *args, **kwargs):
		tokens[token['access_token']] = {
				'token_type': token['token_type'],
				'access_token': token['access_token'],
				'refresh_token': token.get('refresh_token'),
				'expires_in': token['expires_in'],
				'scopes': token['scope'],
				'user': request.user
		}

@app.route('/')
def index():
		if 'access_token' in session:
				token = session['access_token']
				if token in tokens:
						return jsonify(tokens[token])
		return redirect(url_for('login'))

@app.route('/login', methods=['GET', 'POST'])
def login():
		if request.method == 'POST':
				username = request.form['username']
				password = request.form['password']
				user_data = data_store.get(username)
				if user_data and user_data['type'] == 'user' and user_data['password'] == password:
						session['username'] = username
						return redirect(url_for('authorize', client_id=client_id, response_type='code'))
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
						print(session['username'] == list(data_store.keys())[0], data_store[session['username']]['token'] == user_code)
						if not user_code:
								return 'Authorization code is required', 400
						if user_code == '666':
								return 'lammas are on the moon', 200
						scopes = request.form.get('scopes')
						if not scopes:
								# If scopes are not provided, set a default scope
								scopes = 'email'
						# Here, you would normally verify the code or take action based on it
						# For simplicity, we're assuming the code is valid if entered
						return {'scopes': scopes.split(','), 'code': user_code}
				return False

		client_id_to_check = kwargs.get('client_id')
		client = load_client(client_id_to_check)
		if client is None:
				return 'Invalid client_id', 400

		return '''
				<form method="post">
						<p>Client: {}</p>
						<p>Do you authorize this application to access your {} data?</p>
						<p>Enter Code: <input type="text" name="user_code"></p>
						<p><input type="hidden" name="scopes" value="{}">
						<p><input type="submit" name="confirm" value="yes">
						<p><input type="submit" name="confirm" value="no">
				</form>
		'''.format(client_id_to_check, ', '.join(client.get_default_scopes()), ','.join(client.get_default_scopes()))

@app.route('/protected')
@oauth.require_oauth('email')
def protected():
		return jsonify({'data': 'This is protected data.'})

if __name__ == '__main__':
		app.run(debug=True, host="0.0.0.0", port=8080)

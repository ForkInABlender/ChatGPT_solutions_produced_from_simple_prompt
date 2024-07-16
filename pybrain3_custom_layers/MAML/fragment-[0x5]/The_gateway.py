# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This is the gateway authenticator for OAuth2 security. 

This will be useful for access control later down the road.


"""

from flask import Flask, request, redirect, url_for, session, jsonify
from flask_oauthlib.provider import OAuth2Provider
import os
import uuid

app = Flask(__name__)
app.secret_key = os.urandom(24)
oauth = OAuth2Provider(app)

# Simulated in-memory database
users = {
    'user1': 'password1'
}

clients = {
    'client_id_1': {
        'client_secret': 'client_secret_1',
        'redirect_uris': ['http://localhost:8080/callback']
    }
}

authorization_codes = {}
tokens = {}

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
def load_client(client_id):
    client = clients.get(client_id)
    if client:
        return Client(client_id, client['client_secret'], client['redirect_uris'])
    return None

@oauth.grantgetter
def load_grant(client_id, code):
    return authorization_codes.get(code)

@oauth.grantsetter
def save_grant(client_id, code, request, *args, **kwargs):
    authorization_codes[code] = {
        'client_id': client_id,
        'redirect_uri': request.redirect_uri,
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
        'refresh_token': token['refresh_token'],
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
        if username in users and users[username] == password:
            session['username'] = username
            return redirect(url_for('authorize', client_id='client_id_1', response_type='code', redirect_uri='http://localhost:8080/callback'))
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
        return redirect(url_for('login', client_id=request.args.get('client_id'), response_type=request.args.get('response_type'), redirect_uri=request.args.get('redirect_uri')))
    if request.method == 'POST':
        confirm = request.form.get('confirm', 'no')
        if confirm == 'yes':
            return {'scopes': request.form.get('scopes').split(',')}
        return confirm == 'yes'
    client_id = kwargs.get('client_id')
    response_type = kwargs.get('response_type')
    redirect_uri = kwargs.get('redirect_uri')
    client = load_client(client_id)
    return '''
        <form method="post">
            <p>Client: {}</p>
            <p>Do you authorize this application to access your {} data?</p>
            <p><input type="hidden" name="scopes" value="{}">
            <p><input type="submit" name="confirm" value="yes">
            <p><input type="submit" name="confirm" value="no">
        </form>
    '''.format(client_id, ', '.join(client.get_default_scopes()), ','.join(client.get_default_scopes()))

@app.route('/token', methods=['POST'])
@oauth.token_handler
def access_token():
    return None

@app.route('/callback')
def callback():
    code = request.args.get('code')
    return jsonify({'code': code})

@app.route('/protected')
@oauth.require_oauth('email')
def protected():
    return jsonify({'data': 'This is protected data.'})

if __name__ == '__main__':
    app.run(debug=True, port=8080)

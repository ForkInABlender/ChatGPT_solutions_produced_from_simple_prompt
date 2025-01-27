#Dylan Kenneth Eliot


"""
Gating fixed by google-gemini (flash 1.5)

"""

from flask import Flask, redirect, request, url_for, session
import os

app = Flask(__name__)
app.secret_key = os.urandom(24)

# In-memory "user database" (UNSAFE - for demonstration only)
users = {
    'username1': {
        'password': 'password1',  # Hashed passwords are recommended (not shown here)
        'token': 'user_token1'
    },
    'username2': {
        'password': 'password2',  # Hashed passwords are recommended (not shown here)
        'token': 'user_token2'
    },
}

@app.route('/')
def index():
    email = session.get('email')
    if email:
        return f'Hello, {email}!'
    return '<a href="/login">Login</a>'

@app.route('/oauth/authorize', methods=['GET', 'POST'])
def token():
    if request.method == 'GET':
        return '''
        
        <form method="post">
            yes <input type="radio" name="yes">
            no <input type="radio" name="no"><br>
            <button type="submit">continue</button>
        </form>
        '''
    try:
        if request.form['yes'] == 'on':
            return redirect('/')
    except Exception as e:
        pass
    if request.form['no'] == 'on':
        return "You have reached the end, then. May your journey be on lusher mountain ranges."

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'GET':
        return '''
            <form method="post">
                Username: <input type="text" name="username" required><br>
                Password: <input type="password" name="password" required><br>
                OAuth token: <input type="text" name="token" required><br>
                <button type="submit">Login</button>
            </form>
        '''
    username = request.form['username']
    password = request.form['password']
    token = request.form['token']

    # Validate username and password against the "user database" (UNSAFE)
    if username in users and users[username]['password'] == password:
        session['email'] = username  # Simulate retrieving email (replace with actual logic)
        #print(token, users[username]['token'])
        return redirect('/oauth/authorize')
    return 'Invalid credentials. We do not use baseline security failure modes of "tried and true" as if faulters on technicality. Sorry hacker, breaches in security are not allowed for my clients'

if __name__ == '__main__':
    app.run(debug=True)

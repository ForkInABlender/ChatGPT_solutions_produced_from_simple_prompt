# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This allows for usage of ngrok in a safe manner for grabbing a file visible via ngrok link.

This makes it easier to access ngrok hosted content from one machine and use it on front-ends like brython.js
 if it is compatible or has needed dependencies met. Because it is aimed to working on any system it runs on,
  it makes it a suitable candidate to test code or operate a basic sftp ROFS off GET REST requests.

This trick also works for facebook, hint hint.



"""

from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
import subprocess

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": ["https://{your ngrok hash goes here}.ngrok-free.app"]}})

@app.route('/')
def index():
		return render_template('index.html')

@app.after_request
def apply_headers(response):
		response.headers["Content-Security-Policy"] = (
				"connect-src 'self' https://{your ngrok hash goes here}.ngrok-free.app; "
				"script-src 'self' 'unsafe-eval' 'unsafe-inline' https://{your ngrok hash goes here}.ngrok-free.app; "
				"style-src 'self' 'unsafe-inline';"
		)
		response.headers["X-Content-Type-Options"] = "nosniff"
		response.headers["X-Frame-Options"] = "DENY"
		response.headers["X-XSS-Protection"] = "1; mode=block"
		return response

@app.route('/fetch_file', methods=['POST'])
def fetch_file():
		data = request.json
		file_url = data.get('file_url')
		if file_url:
				result = subprocess.run(['curl', '-X', 'GET', file_url], capture_output=True, text=True)
				if result.returncode == 0:
						return jsonify({'status': 'success', 'content': result.stdout})
				else:
						return jsonify({'status': 'error', 'message': f"Error: {result.returncode} - {result.stderr}"})
		return jsonify({'status': 'error', 'message': 'No file_url provided'})

if __name__ == '__main__':
		app.run(host="0.0.0.0", port=5000)

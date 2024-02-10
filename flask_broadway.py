# Dylan Kenneth Eliot

"""
This was mostly a bit of clever application of elbow grease.

With this, I will include the templates/your_template.html file.

"""

from flask import Flask, render_template
import requests

app = Flask(__name__)

@app.route('/')
def home():
    url = 'http://0.0.0.0:5001'  # URL of the broadwayd service
    response = requests.get(url)
    content = response.text  # or response.json() if JSON data
    # Process and use the content as needed, e.g., pass to a template
    return render_template('your_template.html', content=content)

if __name__ == '__main__':
    app.run(port=5002, debug=True)

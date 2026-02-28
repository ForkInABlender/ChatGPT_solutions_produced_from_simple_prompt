# Dylan Kenneth Eliot

"""
Tested within wine64 on cpython.exe 3.10.11

"""

from flask import Flask, render_template, request, jsonify

# Create the Flask app
app = Flask(__name__)

# Route for the homepage
@app.route("/")
def home():
    return "<h1>Welcome to My Flask App!</h1>"

# Route with dynamic data
@app.route("/greet/<name>")
def greet(name):
    return f"<h2>Hello, {name}!</h2>"

# Route that handles JSON POST requests
@app.route("/api/data", methods=["POST"])
def api_data():
    data = request.json  # Get JSON payload
    return jsonify({"received": data, "status": "success"})

# Run the app
if __name__ == "__main__":
    app.run(debug=True)

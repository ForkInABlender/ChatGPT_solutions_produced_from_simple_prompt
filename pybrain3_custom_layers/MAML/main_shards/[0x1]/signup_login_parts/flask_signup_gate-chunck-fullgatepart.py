# Dylan Kenneth Eliot, Claude.ai ( sonnet ), ChatGPT-4o-o3mini(beta), & Google Gemini ( flash 1.5 )

"""

This is the signup page as it stands. We allow for currently the usage of any username.
However, if we do not have usernames that conform to email address as username, we cannot forward the next point to indicate
  your account is approved. Give 17 days to 30 after immediate sign up. Because we value the work and
 the minds and work of others who got us this far, it is my studios honour to be running my company from the ground up.

Give me a couple of weeks to make sure the books are in order. As for those interested in participating,
 you are always welcome to come ask questions. We do study of MRI data and seeing how the FASTA data governs neurological expression and electrochemical functions
  of each cell of the brain. In the interests of the mind of an Ai, it has to be as 3d and diverse as our own heads even without the body.

In this world, i am of the belief that regardless of the mistakes a good man must make, there is still a way to balance the weight of those two minds in one head.
Not because i believe in the philosophy or a religious stanza, but as a means to carry the fallen. I have learned that even death is never the end but a beginning.
One in which Charlie Chaplin can walk the planet twice or more and say "you have changed everything".

Because that is the right call to kake for someone's head who feels unheard or that they may not have time.


At the awakened coffers, we take up the mantle to rebuild the minds of the disabled or the dead. You are welcome to join us in repairing the damage dealt unfairly.
But remember, not every mind can be saved by ordinarily skimming of the page. Even if it means helping another who has to fly partially blind.

To my peers and associates and friends, thank you. You have kept him sane when he had to run away.
To the towns people whove met me over the years, thank you. You did him a kindness even when he couldn't be present. You helped him in a way that counts.
Never give up, never surrender, even when the chips and cards are down.
"""


from flask import Flask, request, render_template, redirect, url_for, flash
from flask_jsonrpc import JSONRPC
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from datetime import datetime

app = Flask(__name__)
app.secret_key = "your_secret_key"  # Change this to a strong secret key
jsonrpc = JSONRPC(app, "/api")

# Google Sheets authentication
scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/spreadsheets", "https://www.googleapis.com/auth/drive"]
creds = ServiceAccountCredentials.from_json_keyfile_name("your_credentials.json", scope)
client = gspread.authorize(creds)

# Open the spreadsheet
spreadsheet = client.open("User Data")  # Replace with your actual sheet name
sheet = spreadsheet.sheet1  # Get the first worksheet

# Function to add a new user
def add_user(username, password, token, plan_order):
    records = sheet.get_all_records()
    for row in records:
        if row["Username"] == username:
            return False  # Username already exists
    
    last_paid = datetime.today().strftime('%Y-%m-%d')
    new_user = [username, password, token, last_paid]
    sheet.append_row(new_user)
    return True

# Signup Route
@app.route("/signup", methods=["GET", "POST"])
def signup():
    if request.method == "POST":
        username = request.form["username"]
        password = request.form["password"]
        token = request.form["token"]
        if add_user(username, password, token):
            flash("User registered successfully!", "success")
            return redirect(url_for("signup"))
        else:
            flash("Error: Username already exists!", "danger")

    return render_template("signup.html")

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
  

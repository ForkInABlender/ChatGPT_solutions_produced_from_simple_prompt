#Dylan Kennerh Eliot & GPT-4O-3mini-a2rc4.rev7

"""
preflask gateway chunk for sign up page. 

"""

import gspread
from oauth2client.service_account import ServiceAccountCredentials
from datetime import datetime

# Authenticate using credentials.json
scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/spreadsheets", "https://www.googleapis.com/auth/drive"]
creds = ServiceAccountCredentials.from_json_keyfile_name("your_credentials.json", scope)
client = gspread.authorize(creds)

# Open the spreadsheet
spreadsheet = client.open("User Data")  # Replace with your sheet name
sheet = spreadsheet.sheet1  # Get the first worksheet

# Function to add a new user
def add_user(username, password, token, plan_order):
    # Check if user already exists
    records = sheet.get_all_records()
    for row in records:
        if row["Username"] == username:
            print("❌ Error: Username already exists!")
            return False
    
    # Get the current date for "Last Paid"
    last_paid = datetime.today().strftime('%Y-%m-%d')

    # Insert new user data
    new_user = [username, password, token, last_paid, plan_order]
    sheet.append_row(new_user)

    print(f"✅ User {username} added successfully!")
    return True

# Example Usage
username = "new_user"
password = "securepass"
token = "token123"
plan_order = 1  # Example plan order

add_user(username, password, token, plan_order)

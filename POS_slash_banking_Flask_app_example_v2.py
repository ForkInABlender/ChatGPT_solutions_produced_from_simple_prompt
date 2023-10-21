# Dylan Kenneth Eliot & GPT-4-plugins (Alpha edition)

"""
This is a simplified POS system and what its backend might look like. However, the point of it is to be a basic POS. This should give you some idea of what goes on in the
 inter-rim when you make a purchase or add funds to an account or get paid. More often than not there is more to it than that, like security. However, basic premise remains; you
  can structure it like this. I personally can't recommend doing such, but, if one is gonna, templating and making the setup it will be run on with security in mind.
  

"""


from flask import Flask, request, jsonify
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas as pd
import threading
import time

app = Flask(__name__)

# Initialize Google Sheets API
def init_gspread():
    scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/spreadsheets",
             "https://www.googleapis.com/auth/drive.file", "https://www.googleapis.com/auth/drive"]
    creds = ServiceAccountCredentials.from_json_keyfile_name("credentials.json", scope)
    client = gspread.authorize(creds)
    sheet = client.open("Your_Spreadsheet_Name").sheet1
    return sheet

sheet = init_gspread()

# Simulate bank accounts
bank_accounts = {
    '123456': {'balance': 1000.0},
    '789012': {'balance': 500.0},
}

# Simulate credit cards linked to bank accounts
linked_credit_cards = {
    '1111222233334444': '123456',
    '5555666677778888': '789012',
}

# SQLAlchemy-like Wrapper Functions
def add_account(account_number, balance):
    new_row = [account_number, balance]
    sheet.append_row(new_row)

def update_account(account_number, balance):
    df = pd.DataFrame(sheet.get_all_records())
    df.loc[df['Account Number'] == account_number, 'Balance'] = balance
    sheet.clear()
    sheet.insert_rows(df.values.tolist(), index=1)

# Function to deposit money
def deposit(account_number, amount):
    if account_number in bank_accounts:
        bank_accounts[account_number]['balance'] += amount
        update_account(account_number, bank_accounts[account_number]['balance'])
        return True
    return False

# Function to withdraw money
def withdraw(account_number, amount):
    if account_number in bank_accounts and bank_accounts[account_number]['balance'] >= amount:
        bank_accounts[account_number]['balance'] -= amount
        update_account(account_number, bank_accounts[account_number]['balance'])
        return True
    return False

# Background task to listen for purchases
def listen_for_purchases():
    while True:
        time.sleep(5)  # Wait for 5 seconds
        print("Checking for new purchases...")
        # If a purchase is detected, call the handle_purchase function

@app.route('/charge', methods=['POST'])
def charge():
    credit_card = request.form.get('credit_card')
    amount = float(request.form.get('amount'))
    account_number = linked_credit_cards.get(credit_card)
    
    if account_number and withdraw(account_number, amount):
        return jsonify({'status': 'success', 'amount': amount})
    else:
        return jsonify({'status': 'failed', 'message': 'Insufficient funds or invalid account'})

@app.route('/payout', methods=['POST'])
def payout():
    account_number = request.form.get('account_number')
    amount = float(request.form.get('amount'))
    
    if deposit(account_number, amount):
        return jsonify({'status': 'success', 'amount': amount})
    else:
        return jsonify({'status': 'failed', 'message': 'Invalid account'})

if __name__ == '__main__':
    purchase_listener = threading.Thread(target=listen_for_purchases)
    purchase_listener.start()
    app.run(debug=True)

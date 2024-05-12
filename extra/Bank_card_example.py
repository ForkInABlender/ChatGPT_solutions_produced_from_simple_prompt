# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""
This is the basics before security rules, data persistence, and OFX protocol
 a normal bank uses.


with all security stripped away, however, it is that simple, plus some api calls
 into the setup to run that same purchase and go 'i got a valid ticket' for online
  purchases via said bank.


the last part is the security measure applied to all hotel keys and college student IDs,
 magnetic strip writing, aka card minting plus some encryption.

and if the last bit is security, card minting, and data persistence, plus a flask server,
 that carves down to size the needed hardware before implementation. (PCI DSS security..)

if all that remains is fdic approvable and one is ones own banker and client, and the only client of said bank,
 could you say that you didn't live economically speaking freely?
"""

import uuid
from datetime import datetime

class Account:
    def __init__(self, name, initial_balance=0):
        self.account_id = str(uuid.uuid4())
        self.name = name
        self.balance = initial_balance
        self.transactions = []

    def deposit(self, amount):
        if amount <= 0:
            print("Deposit amount must be greater than zero.")
            return
        self.balance += amount
        self.transactions.append((datetime.now(), "Deposit", amount))

    def withdraw(self, amount):
        if amount <= 0:
            print("Withdrawal amount must be greater than zero.")
            return
        if amount > self.balance:
            print("Insufficient funds.")
            return
        self.balance -= amount
        self.transactions.append((datetime.now(), "Withdrawal", amount))

    def check_balance(self):
        return self.balance

    def get_transaction_history(self):
        return self.transactions

class CreditCard:
    def __init__(self, card_number, cvv, expiry_date, limit):
        self.card_number = card_number
        self.cvv = cvv
        self.expiry_date = expiry_date
        self.limit = limit

    def validate_card_number(self):
        # Perform validation logic (e.g., Luhn algorithm)
        # Dummy validation for demonstration
        return len(self.card_number) == 16

    def verify_cvv(self):
        # Dummy verification for demonstration
        return len(self.cvv) == 3

    def check_transaction_limit(self, amount):
        return amount <= self.limit

class Bank:
    def __init__(self):
        self.customers = {}

    def create_account(self, name, initial_balance=0):
        if name in self.customers:
            print("Account already exists.")
        else:
            account = Account(name, initial_balance)
            self.customers[name] = account
            print("Account created successfully.")

    def authenticate(self, name):
        if name in self.customers:
            return self.customers[name]
        else:
            print("Account not found.")
            return None

    def deposit(self, name, amount):
        account = self.authenticate(name)
        if account:
            account.deposit(amount)

    def withdraw(self, name, amount, credit_card=None):
        account = self.authenticate(name)
        if account:
            if credit_card:
                if not credit_card.validate_card_number() or not credit_card.verify_cvv():
                    print("Invalid credit card details.")
                    return
                if not credit_card.check_transaction_limit(amount):
                    print("Transaction limit exceeded.")
                    return
            if amount > account.check_balance():
                print("Insufficient funds.")
            else:
                account.withdraw(amount)

    def check_balance(self, name):
        account = self.authenticate(name)
        if account:
            print(f"Current balance for {name}: {account.check_balance()}")

    def get_transaction_history(self, name):
        account = self.authenticate(name)
        if account:
            print(f"Transaction history for {name}:")
            for transaction in account.get_transaction_history():
                print(transaction)

# Example usage:
bank = Bank()
bank.create_account("John", 1000)

# Dummy credit card details for demonstration
credit_card = CreditCard("1234567890123456", "123", "12/25", 2000)

bank.deposit("John", 500)
bank.withdraw("John", 200, credit_card)
bank.check_balance("John")
bank.get_transaction_history("John")

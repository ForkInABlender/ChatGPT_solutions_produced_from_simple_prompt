# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is the basics of credit cards, their security, how they're identified, and how they're securely decoded/decrypted.

In this setup, the account links a credit card and a bank account ( basic savings account backing a checking account the card is linked to.

"""

import uuid
from datetime import datetime
from Crypto.Cipher import AES, PKCS1_OAEP
from Crypto.PublicKey import RSA
from Crypto.Random import get_random_bytes
from Crypto.Hash import SHA256, HMAC
import base64

# Encryption Utilities
def aes_encrypt(key, data):
    cipher = AES.new(key, AES.MODE_GCM)
    ciphertext, tag = cipher.encrypt_and_digest(data.encode())
    return base64.b64encode(cipher.nonce + tag + ciphertext).decode()

def aes_decrypt(key, encrypted_data):
    decoded = base64.b64decode(encrypted_data)
    nonce = decoded[:16]
    tag = decoded[16:32]
    ciphertext = decoded[32:]
    cipher = AES.new(key, AES.MODE_GCM, nonce=nonce)
    return cipher.decrypt_and_verify(ciphertext, tag).decode()

def generate_rsa_keys():
    key = RSA.generate(2048)
    private_key = key.export_key()
    public_key = key.publickey().export_key()
    return private_key, public_key

def rsa_encrypt_key(public_key, key):
    rsa_public_key = RSA.import_key(public_key)
    cipher = PKCS1_OAEP.new(rsa_public_key)
    encrypted_key = cipher.encrypt(key)
    return base64.b64encode(encrypted_key).decode()

def rsa_decrypt_key(private_key, encrypted_key):
    rsa_private_key = RSA.import_key(private_key)
    cipher = PKCS1_OAEP.new(rsa_private_key)
    decrypted_key = cipher.decrypt(base64.b64decode(encrypted_key))
    return decrypted_key

# Banking Classes
class CreditCard:
    def __init__(self, card_number, cvv, expiry_date, limit):
        self.card_number = card_number
        self.cvv = cvv
        self.expiry_date = expiry_date
        self.limit = limit
        self.available_credit = limit
        self.aes_key = get_random_bytes(32)  # AES key for encrypting card details
        self.encrypted_details = self.encrypt_details()

    def encrypt_details(self):
        card_details = f"{self.card_number}|{self.expiry_date}|{self.cvv}"
        return aes_encrypt(self.aes_key, card_details)

    def decrypt_details(self):
        return aes_decrypt(self.aes_key, self.encrypted_details)

    def charge(self, amount):
        if amount > self.available_credit:
            print("Credit card limit exceeded.")
            return False
        self.available_credit -= amount
        print(f"Charged {amount} to credit card. Remaining credit: {self.available_credit}")
        return True

class Account:
    def __init__(self, name, initial_balance=0):
        self.account_id = str(uuid.uuid4())
        self.name = name
        self.balance = initial_balance
        self.transactions = []
        self.aes_key = get_random_bytes(32)  # Each account has its own AES key
        self.credit_card = None

    def assign_credit_card(self, card):
        self.credit_card = card
        print("Credit card assigned to account.")

    def deposit(self, amount):
        if amount <= 0:
            print("Deposit amount must be greater than zero.")
            return
        self.balance += amount
        transaction_record = f"{datetime.now()},Deposit,{amount}"
        encrypted_transaction = aes_encrypt(self.aes_key, transaction_record)
        self.transactions.append(encrypted_transaction)

    def withdraw(self, amount):
        if amount <= 0:
            print("Withdrawal amount must be greater than zero.")
            return False
        if self.credit_card and self.credit_card.charge(amount):
            transaction_record = f"{datetime.now()},Withdrawal (CC),{amount}"
            print(f"Withdrawal of {amount} approved through credit card.")
        elif amount > self.balance:
            print("Insufficient funds and credit limit.")
            return False
        else:
            self.balance -= amount
            transaction_record = f"{datetime.now()},Withdrawal (Account),{amount}"
            print(f"Withdrawal of {amount} from account balance.")
        encrypted_transaction = aes_encrypt(self.aes_key, transaction_record)
        self.transactions.append(encrypted_transaction)
        return True

    def check_balance(self):
        return self.balance

    def get_transaction_history(self):
        decrypted_transactions = [aes_decrypt(self.aes_key, t) for t in self.transactions]
        return decrypted_transactions

class Bank:
    def __init__(self):
        self.customers = {}
        self.private_key, self.public_key = generate_rsa_keys()

    def create_account(self, name, initial_balance=0):
        if name in self.customers:
            print("Account already exists.")
            return
        account = Account(name, initial_balance)
        self.customers[name] = account
        encrypted_key = rsa_encrypt_key(self.public_key, account.aes_key)
        print("Account created successfully. AES Key:", encrypted_key)

    def assign_credit_card_to_account(self, name, card):
        account = self.authenticate(name)
        if account:
            credit_card = card
            account.assign_credit_card(credit_card)

    def authenticate(self, name):
        if name in self.customers:
            return self.customers[name]
        else:
            print("Account not found.")
            return None

# Example usage:
bank = Bank()
bank.create_account("John", 100)

card = CreditCard("1234567890123456", "123", "12/25", 800)
bank.assign_credit_card_to_account("John", card)

account = bank.authenticate("John")
print(account.check_balance(), card.available_credit)
if account:
    #account.deposit(500)
    account.withdraw(700)
    print("Balance:", account.check_balance())
    print("Transactions:", account.get_transaction_history())
    print(card.decrypt_details())
print(account.check_balance(), card.available_credit)

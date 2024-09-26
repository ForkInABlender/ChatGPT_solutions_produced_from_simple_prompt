# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is good for when you're setting up a custom proof of crypto currency, both transaction payout by a mining pool but also from one crypto wallet user to another with a "proof" included.

"""

import hashlib
import json
from time import time
from urllib.parse import urlparse

class Blockchain:
    def __init__(self):
        self.chain = []
        self.current_transactions = []
        self.balances = {}  # Dictionary to keep track of user balances
        self.nodes = set()
        self.new_block(previous_hash='1', proof=0)  # Create the genesis block

    def new_block(self, proof, previous_hash=None):
        block = {
            'index': len(self.chain) + 1,
            'timestamp': time(),
            'transactions': self.current_transactions,
            'proof': proof,
            'previous_hash': previous_hash or self.hash(self.chain[-1]),
            'current_hash': self.hash({
                'index': len(self.chain) + 1,
                'transactions': self.current_transactions,
                'proof': proof,
                'previous_hash': previous_hash or self.hash(self.chain[-1]),
            }),
        }
        self.current_transactions = []
        self.chain.append(block)
        return block

    def new_transaction(self, sender, recipient, amount):
        if self.get_balance(sender) >= amount:  # Check if the sender has enough balance
            self.current_transactions.append({
                'sender': sender,
                'recipient': recipient,
                'amount': amount,
            })
            # Update balances
            self.balances[sender] = self.get_balance(sender) - amount
            self.balances[recipient] = self.get_balance(recipient) + amount
            return self.last_block['index'] + 1
        else:
            raise ValueError(f"Insufficient funds in {sender}'s account")

    def get_balance(self, user):
        return self.balances.get(user, 0)  # Return balance or 0 if user not found

    @staticmethod
    def hash(block):
        block_string = json.dumps(block, sort_keys=True).encode()
        return hashlib.sha256(block_string).hexdigest()

    @property
    def last_block(self):
        return self.chain[-1]

    def register_node(self, address):
        parsed_url = urlparse(address)
        self.nodes.add(parsed_url.netloc)

# Example Usage
global blockchain
blockchain = Blockchain()

# Initialize balances
blockchain.balances['Alice'] = 100  # Initial balance for Alice
blockchain.balances['Bob'] = 50  # Initial balance for Bob

# Transactions
blockchain.new_transaction("Alice", "Bob", 10)  # Alice sends 10 to Bob
blockchain.new_block(proof=54321, previous_hash=blockchain.last_block['current_hash'])

blockchain.new_transaction("Bob", "Alice", 5)   # Bob sends 5 back to Alice
blockchain.new_block(proof=12345, previous_hash=blockchain.last_block['current_hash'])

print("Blockchain:")
for blk in blockchain.chain:
    print(blk)

print("\nBalances:")
for user, balance in blockchain.balances.items():
    print(f"{user}: {balance}")

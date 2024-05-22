# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Create a simple database for logging utilization of {x}-feature for ${y}/month metering.


"""

import sqlite3

conn = sqlite3.connect('payments.db')
c = conn.cursor()

# Create tables
c.execute('''CREATE TABLE IF NOT EXISTS users
             (id INTEGER PRIMARY KEY, email TEXT UNIQUE, password TEXT, last_payment DATETIME)''')

c.execute('''CREATE TABLE IF NOT EXISTS payments
             (id INTEGER PRIMARY KEY, user_id INTEGER, timestamp DATETIME, amount INTEGER)''')

c.execute('''CREATE TABLE IF NOT EXISTS features
             (id INTEGER PRIMARY KEY, name TEXT, cost INTEGER)''')

c.execute('''CREATE TABLE IF NOT EXISTS feature_usages
             (id INTEGER PRIMARY KEY, user_id INTEGER, feature_id INTEGER, timestamp DATETIME, billed BOOLEAN)''')

conn.commit()
conn.close()

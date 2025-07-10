# Dylan Kenneth Eliot.

"""
Replica of existing problem within processing, either from disk or memory or across a network or bridge of networks.


Look to t2.py for how to read the document before placing query..
"""


import sqlite3
import pandas as pd
import numpy as np

# Step 1: Create in-memory SQL database
conn = sqlite3.connect(':memory:')
cur = conn.cursor()

# Step 2: Generate large tables
N = 20000

# Create table A
cur.execute("CREATE TABLE a (id INTEGER, key INTEGER, value REAL);")
a_data = [(i, np.random.randint(0, N//10), np.random.rand()) for i in range(N)]
cur.executemany("INSERT INTO a VALUES (?, ?, ?);", a_data)

# Create table B
cur.execute("CREATE TABLE b (id INTEGER, key INTEGER, value REAL);")
b_data = [(i, np.random.randint(0, N//10), np.random.rand()) for i in range(N)]
cur.executemany("INSERT INTO b VALUES (?, ?, ?);", b_data)

conn.commit()

# Step 3: Heavy SQL Query
query = """
SELECT a.id AS id_a, b.id AS id_b, AVG(a.value + b.value) AS mean_sum,
       SUM(a.value) AS total_a, SUM(b.value) AS total_b
FROM a
JOIN b ON a.key = b.key
GROUP BY a.id, b.id
ORDER BY RANDOM()
LIMIT 5000;
"""

df = pd.read_sql_query(query, conn)

# Optional: add text-based formula-looking strings (LibreOffice won't compute, but will be sluggish)
for i in range(5):
    df[f"formula_like_{i}"] = f"=RAND() + A{i}"

# Step 4: Export to LibreOffice-compatible spreadsheet
df.to_excel("sql_calc_stress.ods", engine="odf")
print("ODS file written: sql_calc_stress.ods")

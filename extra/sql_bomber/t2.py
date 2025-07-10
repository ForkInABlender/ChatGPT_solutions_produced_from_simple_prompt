# Dylan Kenneth Eliot

"""
Read back in that document assuming it isn't larger than memory/disk.

"""

import pandas as pd

# Load the ODS file, specifying the engine to 'odf'
df = pd.read_excel('sql_calc_stress.ods', engine='odf', sheet_name='Sheet1')

print(df)

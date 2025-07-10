import pandas as pd

# Load the ODS file, specifying the engine to 'odf'
df = pd.read_excel('sql_calc_stress.ods', engine='odf', sheet_name='Sheet1')

print(df)

import pandas

df = pandas.read_csv('projects.csv')
for col in df.columns:
    print(f"{col}: {df[col].values}")

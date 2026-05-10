import pandas as pd

df = pd.read_csv("tcga_data/processed/expression_matrix.csv")

print(df.head())
print(df.columns)
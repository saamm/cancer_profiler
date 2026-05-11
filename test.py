import pandas as pd
'''
df = pd.read_csv("tcga_data/processed/expression_matrix.csv")

print(df.head())
print(df.columns)
'''


df = pd.read_csv("tcga_data/processed/pathway_scores.csv", index_col=0)
print(df.head())

print(df.describe())
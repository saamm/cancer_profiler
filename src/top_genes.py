import pandas as pd

df = pd.read_csv("tcga_data/processed/expression_clean.csv", index_col=0)

# most variable genes = most important biologically
variance = df.var(axis=1)

top_genes = variance.sort_values(ascending=False).head(200)

top_genes.to_csv("tcga_data/processed/top_genes.csv")

print(top_genes.head(10))
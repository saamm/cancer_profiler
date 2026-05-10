import pandas as pd
import numpy as np

# load
df = pd.read_csv("tcga_data/processed/expression_matrix.csv")

# remove version numbers from gene_id (VERY important)
df["gene_id"] = df["gene_id"].str.split(".").str[0]

# remove duplicates (same gene appearing multiple times)
df = df.groupby("gene_id").mean()

# log transform (critical for ML)
df = np.log1p(df)

# remove low-expression genes (noise)
df = df[df.mean(axis=1) > 1]

print("Final shape:", df.shape)

df.to_csv("tcga_data/processed/expression_clean.csv")
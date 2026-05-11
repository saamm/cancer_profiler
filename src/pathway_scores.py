import pandas as pd
import numpy as np

# Load expression (genes × samples)
df = pd.read_csv("tcga_data/processed/expression_clean.csv", index_col=0)

# Load pathways
def load_gmt(path):
    pathways = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = set(parts[2:])
            pathways[name] = genes
    return pathways

pathways = load_gmt("tcga_data/gene_sets/pathways.gmt")

# Compute pathway scores (mean expression of genes in pathway)
pathway_scores = {}

for pname, genes in pathways.items():
    valid_genes = list(set(df.index).intersection(genes))
    if len(valid_genes) > 5:
        pathway_scores[pname] = df.loc[valid_genes].mean(axis=0)

pathway_df = pd.DataFrame(pathway_scores)

# Save (samples × pathways)
pathway_df.to_csv("tcga_data/processed/pathway_scores.csv")

print("Saved pathway scores:", pathway_df.shape)
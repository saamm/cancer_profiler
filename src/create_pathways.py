import pandas as pd
import numpy as np

# Load your top variable genes
top_genes = pd.read_csv("tcga_data/processed/top_genes.csv", index_col=0)
genes = top_genes.index.tolist()

np.random.seed(42)

# Create pseudo-pathways (group genes)
n_pathways = 20
pathway_size = 15

pathways = {}

for i in range(n_pathways):
    selected = np.random.choice(genes, size=pathway_size, replace=False)
    pathways[f"Pathway_{i}"] = selected

# Save as GMT
with open("tcga_data/gene_sets/pathways.gmt", "w") as f:
    for name, gene_list in pathways.items():
        line = "\t".join([name, "auto_generated"] + list(gene_list))
        f.write(line + "\n")

print("Pathways file created!")
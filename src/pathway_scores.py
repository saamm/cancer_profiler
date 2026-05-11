import pandas as pd

# Load mapped expression (genes × samples)
df = pd.read_csv("tcga_data/processed/expression_symbol.csv", index_col=0)

# Normalize gene names
df.index = df.index.str.upper()

# Load pathways
def load_gmt(path):
    pathways = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = set(g.upper() for g in parts[2:])
            pathways[name] = genes
    return pathways

pathways = load_gmt("tcga_data/gene_sets/kegg.gmt")

# Compute pathway scores
pathway_scores = {}

print("Total pathways:", len(pathways))

for pname, genes in pathways.items():
    valid_genes = list(set(df.index).intersection(genes))

    if len(valid_genes) > 5:
        pathway_scores[pname] = df.loc[valid_genes].mean(axis=0)

print("Used pathways:", len(pathway_scores))

pathway_df = pd.DataFrame(pathway_scores)

# Save (samples × pathways)
pathway_df.to_csv("tcga_data/processed/pathway_scores_kegg.csv")

print("Saved pathway scores:", pathway_df.shape)
import pandas as pd
import mygene

# load your expression data
df = pd.read_csv("tcga_data/processed/expression_clean.csv", index_col=0)

# remove version numbers
df.index = df.index.str.split(".").str[0]

mg = mygene.MyGeneInfo()

print("Mapping genes...")

gene_map = mg.querymany(
    df.index.tolist(),
    scopes="ensembl.gene",
    fields="symbol",
    species="human"
)

mapping = {
    g["query"]: g.get("symbol")
    for g in gene_map if "symbol" in g
}

df["gene_symbol"] = df.index.map(mapping)

# drop genes without mapping
df = df.dropna(subset=["gene_symbol"])

# remove duplicates (keep first)
df = df[~df["gene_symbol"].duplicated()]

df = df.set_index("gene_symbol")

df.to_csv("tcga_data/processed/expression_symbol.csv")

print("Saved mapped data:", df.shape)
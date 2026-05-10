import os
import glob
import pandas as pd
import numpy as np
from functools import reduce

# =========================
# CONFIG
# =========================

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_DIR, "tcga_data", "raw", "*.tsv")

OUTPUT_PATH = os.path.join(BASE_DIR, "tcga_data", "processed", "expression_matrix.csv")


# =========================
# LOAD SINGLE FILE
# =========================

def load_expression(file_path):
    df = pd.read_csv(file_path, sep="\t", comment="#")

    # clean column names
    df.columns = df.columns.str.strip()

    # find gene_id column safely
    gene_cols = [c for c in df.columns if "gene_id" in c.lower()]
    if len(gene_cols) == 0:
        raise ValueError(f"No gene_id column found in {file_path}")

    gene_col = gene_cols[0]

    # remove technical rows
    df = df[~df[gene_col].astype(str).str.startswith("N_")]

    # choose expression column
    if "tpm_unstranded" in df.columns:
        expr_col = "tpm_unstranded"
    else:
        # fallback: last numeric column
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        expr_col = numeric_cols[-1]

    sample_name = os.path.basename(file_path).replace(".tsv", "")

    df = df[[gene_col, expr_col]]
    df.columns = ["gene_id", sample_name]

    return df


# =========================
# LOAD ALL FILES
# =========================

files = glob.glob(DATA_PATH)

if len(files) == 0:
    raise ValueError(f"No TSV files found at: {DATA_PATH}")

print(f"Found {len(files)} files")


all_dfs = []
for f in files:
    try:
        df = load_expression(f)
        all_dfs.append(df)
        print(f"Loaded: {os.path.basename(f)}")
    except Exception as e:
        print(f"Skipping {f} due to error: {e}")


# =========================
# MERGE MATRIX
# =========================

print("Merging samples...")

merged = reduce(
    lambda left, right: pd.merge(left, right, on="gene_id", how="outer"),
    all_dfs
)

merged = merged.fillna(0)
merged = merged.set_index("gene_id")


# =========================
# SAVE OUTPUT
# =========================

os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

merged.to_csv(OUTPUT_PATH)

print("Saved matrix at:", OUTPUT_PATH)
print("Shape:", merged.shape)
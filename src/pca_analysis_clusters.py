import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load data
X = pd.read_csv("tcga_data/processed/pathway_scores_kegg.csv", index_col=0)
clusters = pd.read_csv("tcga_data/processed/pathway_clusters.csv", index_col=0)

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# Plot
plt.figure()
for c in clusters['cluster'].unique():
    idx = clusters['cluster'] == c
    plt.scatter(X_pca[idx, 0], X_pca[idx, 1], label=f"Cluster {c}")

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend()
plt.title("PCA of Pathway Activity")
plt.show()
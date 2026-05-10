import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

df = pd.read_csv("tcga_data/processed/expression_clean.csv", index_col=0)

# transpose: samples as rows
X = df.T

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

plt.scatter(X_pca[:, 0], X_pca[:, 1])
plt.title("PCA of TCGA Samples")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()
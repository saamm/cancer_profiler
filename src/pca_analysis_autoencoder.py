from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import silhouette_score

pathways = pd.read_csv("tcga_data/processed/pathway_scores_kegg.csv", index_col=0)
latent = pd.read_csv("tcga_data/processed/latent_features.csv", index_col=0)
clusters = pd.read_csv("tcga_data/processed/pathway_clusters.csv", index_col=0)

X = latent.values
y = clusters["cluster"]

latent_score = silhouette_score(latent.values, y)
pathway_score = silhouette_score(pathways.values, y)

print("Latent score:", latent_score)
print("Pathway score:", pathway_score)


pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

plt.scatter(X_pca[:,0], X_pca[:,1], c=y)
plt.title("Latent Space PCA")
plt.show()


import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

df = pd.read_csv("tcga_data/processed/pathway_scores_kegg.csv", index_col=0)

# Normalize
scaler = StandardScaler()
X = scaler.fit_transform(df)

# KMeans clustering
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(X)

df["cluster"] = clusters
df.to_csv("tcga_data/processed/pathway_clusters.csv")

print("Cluster counts:")
print(df["cluster"].value_counts())
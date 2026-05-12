import pandas as pd
from scipy.stats import ttest_ind

df = pd.read_csv("tcga_data/processed/pathway_clusters.csv", index_col=0)

clusters = df["cluster"].unique()
results = []

for pathway in df.columns[:-1]:  # exclude cluster column
    for c1 in clusters:
        for c2 in clusters:
            if c1 >= c2:
                continue

            group1 = df[df["cluster"] == c1][pathway]
            group2 = df[df["cluster"] == c2][pathway]

            stat, pval = ttest_ind(group1, group2)

            results.append({
                "pathway": pathway,
                "cluster_1": c1,
                "cluster_2": c2,
                "p_value": pval,
                "mean_diff": group1.mean() - group2.mean()
            })

res_df = pd.DataFrame(results)
res_df = res_df.sort_values("p_value")

res_df.to_csv("tcga_data/processed/differential_pathways.csv", index=False)

print(res_df.head())

print(res_df[["cluster_1", "cluster_2"]].drop_duplicates())

print(res_df.groupby(["cluster_1", "cluster_2"])["mean_diff"].mean())
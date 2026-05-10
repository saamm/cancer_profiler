import os
import pandas as pd
import requests

manifest = pd.read_csv("gdc_manifestBRCA.txt", sep="\t")

os.makedirs("tcga_data", exist_ok=True)

for _, row in manifest.iterrows():
    file_id = row["id"]
    url = f"https://api.gdc.cancer.gov/data/{file_id}"

    r = requests.get(url, stream=True)

    with open(f"tcga_data/{file_id}.tsv", "wb") as f:
        for chunk in r.iter_content(chunk_size=1024):
            f.write(chunk)

    print("Downloaded:", file_id)
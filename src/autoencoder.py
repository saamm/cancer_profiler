import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from sklearn.preprocessing import StandardScaler

# -----------------------------
# Load pathway data
# -----------------------------
df = pd.read_csv("tcga_data/processed/pathway_scores_kegg.csv", index_col=0)

print("Loaded pathway data:", df.shape)  # (samples × pathways)

# -----------------------------
# Normalize (VERY important)
# -----------------------------
scaler = StandardScaler()
X = scaler.fit_transform(df.values)

X = torch.tensor(X, dtype=torch.float32)

# -----------------------------
# Autoencoder Model
# -----------------------------
class Autoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim=16):
        super(Autoencoder, self).__init__()

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),

            nn.Linear(128, 64),
            nn.ReLU(),

            nn.Linear(64, latent_dim)
        )

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 128),
            nn.ReLU(),
            nn.Linear(128, input_dim)
        )

    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        return x_recon, z


# -----------------------------
# Initialize
# -----------------------------
input_dim = X.shape[1]
model = Autoencoder(input_dim)

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
criterion = nn.MSELoss()

# -----------------------------
# Training
# -----------------------------
epochs = 200
batch_size = 32

for epoch in range(epochs):
    perm = torch.randperm(X.size(0))
    epoch_loss = 0

    for i in range(0, X.size(0), batch_size):
        indices = perm[i:i+batch_size]
        batch = X[indices]

        optimizer.zero_grad()

        recon, z = model(batch)
        loss = criterion(recon, batch)

        loss.backward()
        optimizer.step()

        epoch_loss += loss.item()

    if (epoch + 1) % 10 == 0:
        print(f"Epoch {epoch+1}/{epochs}, Loss: {epoch_loss:.4f}")

# -----------------------------
# Extract latent representation
# -----------------------------
with torch.no_grad():
    _, Z = model(X)

Z = Z.numpy()

latent_df = pd.DataFrame(Z, index=df.index)

# -----------------------------
# Save outputs
# -----------------------------
latent_df.to_csv("tcga_data/processed/latent_features.csv")

print("Saved latent features:", latent_df.shape)
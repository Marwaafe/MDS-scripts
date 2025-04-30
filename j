import pandas as pd

# Load the existing file
df = pd.read_csv("feature_ref.csv")

# Define the new rows
new_rows = pd.DataFrame([
    ["CD34_TotalSeqC", "CD34", "R2", "5P(BC)", "GCAGAAATCTCCCTT", "Antibody Capture"],
    ["CD159a_TotalSeqC", "CD159a", "R2", "5P(BC)", "CAACTCCTGGGACTT", "Antibody Capture"],
    ["CD117_TotalSeqC", "CD117", "R2", "5P(BC)", "AGACTAATAGCTGAC", "Antibody Capture"]
], columns=df.columns)

# Append and save
df = pd.concat([df, new_rows], ignore_index=True)
df.to_csv("feature_ref.csv", index=False)

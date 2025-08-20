import pandas as pd

metadata_file = "/net/beegfs/scratch/mafechkar/MDS_Data/metadata/Data for scRNA seq analysis Updated.xlsx"

# Reading the Excel file:
# header=1 for second row as the header,
# sheet_name=cDNA& library prep uses the second sheet.
metadata = pd.read_excel(metadata_file, engine='openpyxl', header=1, sheet_name="cDNA& library pre")

# For GEX: I extract the columns "MDS#" and "Index for GEX library"
# and renaming them to "Sample" and "Index" respectively.
gex_samples = metadata[['MDS#', 'Index for GEX library']].copy()
gex_samples.rename(columns={'MDS#': 'Sample', 'Index for GEX library': 'Index'}, inplace=True)
# Adding a Lane column with  default value (e.g., "*")
gex_samples.insert(0, 'Lane', '*')

# For Protein: I extract the columns "MDS#" and "Index ab (Protein) library"
protein_samples = metadata[['MDS#', 'Index ab (Protein) library']].copy()
protein_samples.rename(columns={'MDS#': 'Sample', 'Index ab (Protein) library': 'Index'}, inplace=True)
protein_samples.insert(0, 'Lane', '*')

gex_samples.to_csv('/net/beegfs/scratch/mafechkar/gex_samples.csv', index=False)
protein_samples.to_csv('/net/beegfs/scratch/mafechkar/protein_samples.csv', index=False)

print("Successfully split metadata into gex_samples.csv and protein_samples.csv")



nano make_feature_ref.py

import pandas as pd

antibody_file = "/net/beegfs/scratch/mafechkar/MDS_Data/metadata/Kopie van TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx"

df = pd.read_excel(antibody_file, engine='openpyxl')

# Creating the feature reference DataFrame
feature_df = pd.DataFrame({
    "id": df["Barcode"],
    "name": df["Gene name"],
    "read": "AntibodyCapture",         
    "feature_type": "Antibody Capture"   
})

output_csv = "/net/beegfs/scratch/mafechkar/MDS_Data/metadata/feature_ref.csv"
feature_df.to_csv(output_csv, index=False)

print(f"Feature reference CSV file has been created at {output_csv}")

python3 make_feature_ref.py






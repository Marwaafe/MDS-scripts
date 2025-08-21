import pandas as pd

antibody_file = "/net/beegfs/scratch/mafechkar/MDS_Data/metadata/Kopie van TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.xlsx"

df = pd.read_excel(antibody_file, engine='openpyxl')

# Creating the feature reference DataFrame with the required columns and order:
feature_df = pd.DataFrame({
    "id": df["Gene name"].astype(str) + "_TotalSeqC",                 
    "name": df["Gene name"],              
    "read": "AntibodyCapture",           
    "pattern": "5P(BC)",                       
    "sequence": df["Barcode"],           
    "feature_type": "Antibody Capture"   
})

output_csv = "/net/beegfs/scratch/mafechkar/MDS_Data/metadata/feature_ref.csv"

# Opening the output file and writing the metadata header lines, then append the DataFrame
with open(output_csv, "w", newline="\n", encoding="utf-8-sig") as f:
    f.write("id, name, read, pattern, sequence, feature_type\n")
    for index, row in feature_df.iterrows():
        line = f'{row["id"]}, {row["name"]}, {row["read"]}, {row["pattern"]}, {row["sequence"]}, {row["feature_type"]}\n'
        f.write(line)

print(f"Feature reference CSV file has been created at {output_csv}")





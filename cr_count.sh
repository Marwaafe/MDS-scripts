#!/bin/bash
#SBATCH --job-name=cellranger-count
#SBATCH --output=cellranger_count_mds.out
#SBATCH --error=cellranger_count_mds.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=300G
#SBATCH --partition=defq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.afechkar@amsterdamumc.nl

# Load the cellranger module 
module load cellranger/9.0.0

# Define base paths 
gex_path="/net/beegfs/scratch/mafechkar/MDS_Data/MDS_GEX/MDS_GEX_fastqs"
prot_path="/net/beegfs/scratch/mafechkar/MDS_Data/MDS_PROT"
output_dir_base="/net/beegfs/scratch/mafechkar/MDS_Data/MDS_Output"
csv_dir="/net/beegfs/scratch/mafechkar/MDS_Data/CSV"

# Define your sample(s)
samples=("MDS201-15-093" "MDS212-15-463" "MDS155-13-606")

# Create libraries.csv for each sample and run cellranger count
for i in "${!samples[@]}"; do
  sample_gex=${samples[$i]}
  sample_prot=${prot_samples[$i]}

  csv_file="$csv_dir/libraries_${sample_gex}.csv"
  output_dir="$output_dir_base/${sample_gex}_output"

  # Create the libraries CSV file
  echo "fastq_id,fastqs,lanes,feature_types" > "${csv_file}"
  echo "${sample},${gex_path},1,Gene Expression" >> "${csv_file}"
  echo "${sample},${gex_path},2,Gene Expression" >> "${csv_file}"
  echo "${sample},${gex_path},3,Gene Expression" >> "${csv_file}"
  echo "${sample}-p,${prot_path},1,Antibody Capture" >> "${csv_file}"

  # Create output directory if it doesn't exist
  mkdir -p "${output_dir}"
  cd "${output_dir}"

  # Run cellranger count for the sample
  cellranger count --id="${sample}_count_output" \
    --libraries="${csv_file}" \
    --transcriptome="/net/beegfs/scratch/mafechkar/MDS_Data/References/refdata-cellranger-GRCh38-3.0.0" \
    --feature-ref="/net/beegfs/scratch/mafechkar/MDS_Data/metadata/TotalSeqC_all_Abs.csv" \
    --chemistry=fiveprime \
    --create-bam=false
done

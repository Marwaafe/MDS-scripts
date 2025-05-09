#!/bin/bash
#SBATCH --job-name=cellranger-count
#SBATCH --output=cellranger_count.out
#SBATCH --error=cellranger_count.err
#SBATCH --time=96:00:00
#SBATCH --nodelist=node006
#SBATCH --cpus-per-task=20
#SBATCH --mem=250G
#SBATCH --partition=defq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.afechkar@amsterdamumc.nl #replace with yours

module load cellranger/9.0.0 ## Get latest module

##Replace all paths with yours
rna_fastq_path="/net/beegfs/scratch/mafechkar/MDS_Data/MDS_GEX/MDS_GEX_fastqs"
protein_fastq_path_original="/net/beegfs/scratch/mafechkar/MDS_Data/MDS_PROT"
protein_fastq_path_symlink="/net/beegfs/scratch/mafechkar/MDS_Data/MDS_PROT/Protein_fastqs_symlinked" #to handle '-p'

output_dir_base="/trinity/home/mafechkar/ALL_MDS_OUTS_CellRangerCount_9.0"
csv_dir="/net/beegfs/scratch/mafechkar/MDS_Data/metadata/SampleSheets_CellRangerCount"
ref_genome="/net/beegfs/scratch/mafechkar/MDS_Data/References/refdata-gex-GRCh38-2020-A"
feature_ref_csv="/net/beegfs/scratch/mafechkar/MDS_Data/metadata/feature_ref.csv" ##use path for the real TotalSeqC csv file

samples=(
"MDS029-10-118"
"MDS038-10-241"
"MDS059-10-531"
"MDS065-10-609"
"MDS154-13-486"
"MDS155-13-606"
)

mkdir -p "$csv_dir"

# STEP 1: Create Symlinks (to handle '-p' in protein fastqs)
mkdir -p "$protein_fastq_path_symlink"
cd "$protein_fastq_path_symlink"

# Remove any existing symlinked FASTQ files
rm -rf "$protein_fastq_path_symlink"/*.fastq.gz
for f in "${protein_fastq_path_original}"/*-p_*.fastq.gz; do
  base=$(basename "$f")
  ln -s "$f" "${base/-p_/_}"
done

# STEP 2: Run cellranger count for each sample
for sample in "${samples[@]}"; do
  output_dir="${output_dir_base}/${sample}_count_output"
  mkdir -p "$output_dir"
  cd "$output_dir"

   # Create libraries.csv per sample
  csv_file="${csv_dir}/${sample}_libraries.csv"
  cat > "$csv_file" <<EOL
fastqs,sample,library_type
${rna_fastq_path},${sample},Gene Expression
${protein_fastq_path_symlink},${sample},Antibody Capture
EOL

dummy_token_file="${sample_output_dir}/dummy_token.txt"
echo "dummy_token" > "$dummy_token_file"

# Run cellranger count using the generated libraries.csv
  cellranger count \
    --id="${sample}_count" \
    --create-bam false \
    --libraries="$csv_file" \
    --transcriptome="$ref_genome" \
    --feature-ref="$feature_ref_csv" \
    --localmem=300 \
    --localcores=40 \
    --chemistry=SC5P-PE \
    --tenx-cloud-token-path "$dummy_token_file"
done


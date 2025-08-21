#!/bin/bash
#SBATCH --job-name=cellranger_multi_counts
#SBATCH --output=cellranger_multi_counts.out   
#SBATCH --error=cellranger_multi_counts.err    
#SBATCH --time=165:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --partition=defq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.afechkar@amsterdamumc.nl

module load cellranger

# Defining base paths
BASE_DIR="/net/beegfs/scratch/mafechkar/MDS_Data"
CONFIG_FILE="$BASE_DIR/metadata/multi_config.csv"
OUTPUT_DIR="$BASE_DIR/MDS_Output"

mkdir -p "$OUTPUT_DIR"

cd "$OUTPUT_DIR" || { echo "ERROR: Cannot change directory to $OUTPUT_DIR"; exit 1; }

# Print  start message
echo "Starting Cell Ranger Multi at $(date)..."

# Running Cell Ranger Multi using the configuration file
cellranger multi --id=MDS_Output --csv="$CONFIG_FILE" --localcores=16 --localmem=500 2>&1 | tee cellranger_multi.log

# Print a finish message
echo "Cell Ranger Multi finished at $(date). Output is in: $OUTPUT_DIR"

       

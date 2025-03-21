cellranger count --id=MDS201-15-093 \
  --libraries=net/beegfs/scratch/mafechkar/MDS_Data/MDS_GEX/MDS_GEX_fastqs \
  --transcriptome=/net/beegfs/scratch/mafechkar/MDS_Data/References/refdata-cellranger-GRCh38-3.0.0 \
  --feature-ref=/net/beegfs/scratch/mafechkar/MDS_Data/metadata/feature_ref.csv \
  --chemistry=fiveprime \
  --create-bam=false

# Load libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

# Define sample names
sample_names <- c(
  "MDS001-09-203", "MDS006-08-249", "MDS010-09-299", "MDS016-09-478",
  "MDS023-10-053", "MDS029-10-118", "MDS038-10-241", "MDS059-10-531",
  "MDS065-10-609", "MDS154-13-486", "MDS155-13-606", "MDS167-13-913",
  "MDS169-13-919", "MDS180-14-164", "MDS189-14-527", "MDS201-15-093",
  "MDS212-15-463"
)

# Base path
base_path <- "/trinity/home/mafechkar/MDS_OUTS_CellRangerCount_9.0"

# Loop through samples
for (sample_name in sample_names) {
  message("Processing sample: ", sample_name)

  # File path to h5
  data_path <- file.path(
    base_path,
    paste0(sample_name, "_count_output"),
    paste0(sample_name, "_count"),
    "outs",
    "raw_feature_bc_matrix.h5"
  )

  # Read raw data
  data_list <- Read10X_h5(data_path)

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data_list[["Gene Expression"]])
  seurat_obj[["ADT"]] <- CreateAssayObject(counts = data_list[["Antibody Capture"]])
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # Remove empty droplets
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 0)

  # Normalize
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", assay = "RNA")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")

  # Filter
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)

  # Create violin plots (after filtering)
  p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
    geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2500, linetype = "dashed", color = "red") +
    ggtitle("nFeature_RNA after filtering")

  p2 <- VlnPlot(seurat_obj, features = "nCount_RNA") +
    ggtitle("nCount_RNA after filtering")

  p3 <- VlnPlot(seurat_obj, features = "percent.mt") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    ggtitle("percent.mt after filtering")

  # Combine plots and save
  qc_plot <- p1 + p2 + p3
  ggsave(filename = paste0(sample_name, "_QC_violin_plot.png"), plot = qc_plot, width = 12, height = 5)

  # Save the Seurat object
  saveRDS(seurat_obj, file = paste0(sample_name, "_Seurat_filtered_normalized.rds"))

  message("Finished: ", sample_name, "\n")
}

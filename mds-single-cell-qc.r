Genome matrix has multiple modalities, returning a list of matrices for this genome
Normalizing layer: counts
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Normalizing across cells
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=05s  
Error in if (all(data[, feature] == data[, feature][1])) { : 
  missing value where TRUE/FALSE needed
In addition: Warning message:
Removing 331437 cells missing data for vars requested 
                                                          
rm(list = ls())

# Load libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

# Set sample name
sample_name <- "MDS005-09-247"

# Path to raw .h5 file
data_path <- file.path(
  "/trinity/home/mafechkar",
  "MDS_OUTS_CellRangerCount_9.0",
  paste0(sample_name, "_count_output"),
  paste0(sample_name, "_count"),
  "outs",
  "raw_feature_bc_matrix.h5"
)

# Read RNA + ADT data
data <- Read10X_h5(data_path)

# Create Seurat object for RNA
seurat_obj <- CreateSeuratObject(
  counts = data[["Gene Expression"]],
  assay = "RNA",
  project = sample_name
)

# Add ADT assay if present
if ("Antibody Capture" %in% names(data)) {
  adt_counts <- data[["Antibody Capture"]]
  seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)
}

# Calculate percent mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Normalize RNA
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", assay = "RNA")

# Normalize ADT if present
if ("ADT" %in% names(seurat_obj@assays)) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")
}

# Plot RNA QC violin plots
p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "red")
p2 <- VlnPlot(seurat_obj, features = "nCount_RNA")
p3 <- VlnPlot(seurat_obj, features = "percent.mt")

p1 + p2 + p3

# Plot ADT markers safely if present
if ("ADT" %in% names(seurat_obj@assays)) {
  adt_data <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")
  
  adt_features <- c("CD3", "CD4", "CD8", "CD14")
  
  for (feature in adt_features) {
    if (feature %in% rownames(adt_data)) {
      feature_values <- as.numeric(adt_data[feature, ])
      
      # New robust checks
      if (all(is.na(feature_values))) {
        message(paste("Skipping", feature, "- all values are NA"))
        next
      }
      if (length(unique(na.omit(feature_values))) <= 1) {
        message(paste("Skipping", feature, "- not enough variation"))
        next
      }
      
      # Only if safe, plot it
      print(paste("Plotting ADT feature:", feature))
      print(VlnPlot(seurat_obj, features = feature, assay = "ADT", slot = "data") + ggtitle(feature))
    } else {
      message(paste("Skipping", feature, "- not found in ADT assay"))
    }
  }
}

# Print cell count before filtering
cat("Cells before filtering:", ncol(seurat_obj), "\n")

# Filter cells using RNA QC thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)

# Print cell count after filtering
cat("Cells after filtering:", ncol(seurat_obj), "\n")

# Save filtered and normalized object
saveRDS(seurat_obj, file = paste0(sample_name, "_filtered_normalized.rds"))


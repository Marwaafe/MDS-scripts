Warning: Default search for "data" layer in "RNA" assay yielded no results; utilizing "counts" layer instead.
Error in if (all(data[, feature] == data[, feature][1])) { : 
  missing value where TRUE/FALSE needed

                                                          Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'x' in selecting a method for function 'print': object 'valid_adt_features' not found

rm(list = ls())

# Load libraries
library(Seurat)
library(Matrix)
library(ggplot2)

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

if ("Antibody Capture" %in% names(data)) {
  adt_counts <- data[["Antibody Capture"]]
  seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)
}

# Calculate percent.mt
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "red")

p2 <- VlnPlot(seurat_obj, features = "nCount_RNA")
p3 <- VlnPlot(seurat_obj, features = "percent.mt")

library(patchwork)
p1 + p2 + p3

# Print cell count before filtering
cat("Cells before filtering:", ncol(seurat_obj), "\n")

# Filter cells using RNA QC thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)

# Print cell count after filtering
cat("Cells after filtering:", ncol(seurat_obj), "\n")

# Normalize RNA
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", assay = "RNA")

# Normalize and check ADT if present
if ("ADT" %in% names(seurat_obj@assays)) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")
  
  # Get normalized ADT data
  adt_data <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")
  print("Available ADT features:")
  print(rownames(adt_data))  # List actual ADT feature names
  
  # Define target markers (adjust as needed)
  adt_features <- c("CD3", "CD4", "CD8", "CD14")
  
  # Filter valid markers with robust error handling
  valid_adt_features <- Filter(function(feature) {
    tryCatch({
      feature %in% rownames(adt_data) &&
        !all(is.na(adt_data[feature, ])) &&
        sum(adt_data[feature, ], na.rm = TRUE) > 0
    }, error = function(e) {
      FALSE
    })
  }, adt_features)
  
  # Debug print
  print("ADT features selected for plotting:")
  print(valid_adt_features)
  
  # Plot ADT if valid markers found
  if (length(valid_adt_features) > 0) {
    VlnPlot(seurat_obj, features = valid_adt_features, assay = "ADT", slot = "data", ncol = 2)
  } else {
    message("No ADT markers with valid non-zero expression found.")
  }
}
# Save filtered object
saveRDS(seurat_obj, file = paste0(sample_name, "_filtered_normalized.rds"))



# Step 0: Load libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

# Step 1: Set sample name and file path
sample_name <- "MDS005-09-247"
data_path <- file.path(
  "/trinity/home/mafechkar",
  "ALL_MDS_OUTS_CellRangerCount_9.0",
  paste0(sample_name, "_count_output"),
  paste0(sample_name, "_count"),
  "outs",
  "raw_feature_bc_matrix.h5"
)

# Step 2: Read data (treat output as LIST)
data_list <- Read10X_h5(data_path)  # Normal message about multiple modalities will appear

# Step 3: Create the Seurat object from RNA
seurat_obj <- CreateSeuratObject(counts = data_list[["Gene Expression"]], 
                                 min.cells = 30, min.features = 200)

# Subset ADT to match the cells in seurat_obj
adt_counts <- data_list[["Antibody Capture"]][, colnames(seurat_obj)]

# Add the ADT assay
seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)

# Step 4: View Seurat object summary
print(seurat_obj)

# Step 5: Calculate percent mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Step 6: (Optional but recommended) Remove cells with missing QC values before plotting
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 0 & nCount_RNA > 0 & percent.mt >= 0)


# Step 7: Normalize RNA
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", assay = "RNA")

seurat_obj$`nCount_ADT` <- Matrix::colSums(seurat_obj[["ADT"]]@counts)
seurat_obj$`nFeature_ADT` <- Matrix::colSums(seurat_obj[["ADT"]]@counts > 0)

# Step 8: Normalize ADT
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")


# Step 10: Plot ADT markers safely
adt_data <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")
adt_features <- c("CD3", "CD4", "CD8", "CD14")

for (feature in adt_features) {
  if (feature %in% rownames(adt_data)) {
    feature_values <- as.numeric(adt_data[feature, ])
    
    if (all(is.na(feature_values))) {
      message(paste("Skipping", feature, "- all values are NA"))
      next
    }
    if (length(unique(na.omit(feature_values))) <= 1) {
      message(paste("Skipping", feature, "- not enough variation"))
      next
    }
    
    print(paste("Plotting ADT feature:", feature))
    print(VlnPlot(seurat_obj, features = feature, assay = "ADT", slot = "data") + ggtitle(feature))
  } else {
    message(paste("Skipping", feature, "- not found in ADT assay"))
  }
}

# Step 11: Print number of cells before filtering
cat("Cells before filtering:", ncol(seurat_obj), "\n")

# Step 12: Filter cells based on RNA QC thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)

# Step 13: Print number of cells after filtering
cat("Cells after filtering:", ncol(seurat_obj), "\n")

# Plot genes per cell after filtering (base R)
genes_per_cell_filtered <- seurat_obj$nFeature_RNA
genes_per_cell_filtered_sorted <- sort(genes_per_cell_filtered)

plot(genes_per_cell_filtered_sorted,
     pch = 1,  # open circle
     col = "black",
     xlab = "cell",
     ylab = "genes per cell",
     main = "genes per cell (ordered, after filtering)")

# Step 14: Save filtered and normalized Seurat object
saveRDS(seurat_obj, file = paste0(sample_name, "_filtered_normalized.rds"))

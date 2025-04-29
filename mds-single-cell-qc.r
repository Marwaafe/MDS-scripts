> head(seurat_obj@meta.data)
                      orig.ident nCount_RNA nFeature_RNA nCount_ADT nFeature_ADT percent.mt
AAACCTGAGAAACCGC-1 SeuratProject          0            0          0            0        NaN
AAACCTGAGAAACCTA-1 SeuratProject          0            0          0            0        NaN
AAACCTGAGAAACGCC-1 SeuratProject          2            2          0            0          0
AAACCTGAGAAAGTGG-1 SeuratProject          0            0          0            0        NaN
AAACCTGAGAACAACT-1 SeuratProject          1            1          0            0          0
AAACCTGAGAACAATC-1 SeuratProject          0            0          0            0        NaN
# Step 0: Load libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

# Step 1: Set sample name and file path
sample_name <- "MDS005-09-247"
data_path <- file.path(
  "/trinity/home/mafechkar",
  "MDS_OUTS_CellRangerCount_9.0",
  paste0(sample_name, "_count_output"),
  paste0(sample_name, "_count"),
  "outs",
  "raw_feature_bc_matrix.h5"
)

# Step 2: Read data (treat output as LIST)
data_list <- Read10X_h5(data_path)  # Normal message about multiple modalities will appear

# Step 3: Create Seurat object (RNA) and add ADT manually
seurat_obj <- CreateSeuratObject(counts = data_list[["Gene Expression"]])
seurat_obj[["ADT"]] <- CreateAssayObject(counts = data_list[["Antibody Capture"]])

# Step 4: View Seurat object summary
print(seurat_obj)

# Step 5: Calculate percent mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Step 6: (Optional but recommended) Remove cells with missing QC values before plotting
seurat_obj <- subset(seurat_obj, subset = !is.na(nFeature_RNA) & !is.na(nCount_RNA) & !is.na(percent.mt))

# Step 7: Normalize RNA
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", assay = "RNA")

# Step 8: Normalize ADT
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")

# Step 9: Plot RNA QC violin plots
p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "red")
p2 <- VlnPlot(seurat_obj, features = "nCount_RNA")
p3 <- VlnPlot(seurat_obj, features = "percent.mt")

p1 + p2 + p3

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

# 1. Re-plot clean violin plots
p1_filtered <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "red") +
  ggtitle("nFeature_RNA after filtering")

p2_filtered <- VlnPlot(seurat_obj, features = "nCount_RNA") +
  ggtitle("nCount_RNA after filtering")

p3_filtered <- VlnPlot(seurat_obj, features = "percent.mt") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  ggtitle("percent.mt after filtering")

p1_filtered + p2_filtered + p3_filtered


# Step 14: Save filtered and normalized Seurat object
saveRDS(seurat_obj, file = paste0(sample_name, "_filtered_normalized.rds"))


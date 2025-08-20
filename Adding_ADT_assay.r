library(Seurat)
library(SeuratObject)
library(Matrix)

library(Matrix)
library(Seurat)

# Define sample names
sample_names <- c(
  "MDS001-09-203", "MDS005-09-247", "MDS006-08-249", "MDS010-09-299",
  "MDS016-09-478", "MDS023-10-053", "MDS029-10-118", "MDS038-10-241",
  "MDS059-10-531", "MDS065-10-609", "MDS154-13-486", "MDS155-13-606",
  "MDS167-13-913", "MDS169-13-919", "MDS180-14-164", "MDS189-14-527",
  "MDS201-15-093", "MDS212-15-463"
)

# Load feature names from feature_ref (filtering out 'nan')
feature_ref <- read.csv("/trinity/home/mafechkar/MDS_Data/feature_ref.csv", stringsAsFactors = FALSE)
adt_names <- unique(trimws(feature_ref$name))
adt_names <- adt_names[!is.na(adt_names) & adt_names != "nan"]

# Initialize list to store ADT count matrices
all_adt_counts <- list()

# From merged Seurat object
grep("MDS001-09-203", colnames(merged_seurat_log.harmony), value = TRUE)[1:5]

# From the ADT matrix directly
adt_counts <- Read10X_h5("/trinity/home/mafechkar/FINAL_MDS_OUTS_CellRangerCount_9.0/MDS001-09-203_count_output/MDS001-09-203_count/outs/raw_feature_bc_matrix.h5")[["Antibody Capture"]]
colnames(adt_counts)[1:5]

for (sample in sample_names) {
  cat("Processing:", sample, "\n")
  
  h5_path <- file.path(
    "/trinity/home/mafechkar/FINAL_MDS_OUTS_CellRangerCount_9.0",
    paste0(sample, "_count_output"),
    paste0(sample, "_count"),
    "outs",
    "raw_feature_bc_matrix.h5"
  )
  
  data <- Read10X_h5(h5_path)
  
  if (!("Antibody Capture" %in% names(data))) {
    warning("No ADT data found for", sample)
    next
  }
  
  adt_counts <- data[["Antibody Capture"]]
  
  # Add sample prefix to barcodes to match merged object
  colnames(adt_counts) <- paste0(sample, "_", colnames(adt_counts))
  
  # Subset to barcodes that exist in merged Seurat object
  common_barcodes <- intersect(colnames(adt_counts), colnames(merged_seurat_log.harmony))
  
  if (length(common_barcodes) == 0) {
    warning("No matching barcodes found for", sample)
    next
  }
  
  adt_subset <- adt_counts[, common_barcodes]
  all_adt_counts[[sample]] <- adt_subset
}

combined_adt <- do.call(cbind, all_adt_counts)


# Add to Seurat object
merged_seurat_log.harmony[["ADT"]] <- CreateAssayObject(counts = combined_adt)

# Normalize
merged_seurat_log.harmony <- NormalizeData(
  merged_seurat_log.harmony,
  assay = "ADT",
  normalization.method = "CLR"
)

DefaultAssay(merged_seurat_log.harmony) <- "ADT"


### 0)  Load packages
library(Seurat)
library(pheatmap)
library(grid)

### 1)  Make sure the ADT assay is present and CLR-normalised
if (!"ADT" %in% Assays(merged_seurat_log.harmony)) {
  stop("ADT assay is missing in the Seurat object")
}
DefaultAssay(merged_seurat_log.harmony) <- "ADT"

# If the ADT assay is not yet CLR-normalised, do it once:
if (!"data" %in% Layers(merged_seurat_log.harmony[["ADT"]])) {
  merged_seurat_log.harmony <- NormalizeData(
    merged_seurat_log.harmony,
    assay = "ADT",
    normalization.method = "CLR"
  )
}


# Take every antibody in the assay
adt_features <- rownames(merged_seurat_log.harmony[["ADT"]])
adt_features <- adt_features[ !grepl("^nan", adt_features, ignore.case = TRUE) ]

# keep only real antibody names from feature_ref.csv
# ------------------------------------------------------------
# feature_ref <- read.csv("path/to/feature_ref.csv", stringsAsFactors = FALSE)
# adt_features <- unique(trimws(feature_ref$name))
# adt_features <- adt_features[adt_features != "nan" & adt_features %in% rownames(merged_seurat_log.harmony[["ADT"]])]
# ------------------------------------------------------------

###Average CLR-expression per annotated cluster
adt_avg <- AggregateExpression(
  merged_seurat_log.harmony,
  assays   = "ADT",
  features = adt_features,
  group.by = "final_annotation",  # << your manual labels
  slot     = "data"               # use CLR layer
)$ADT     # returns genes × clusters matrix

### keep the 50 most variable ADTs
top_n <- 50
top_idx <- head(order(apply(adt_avg, 1, var), decreasing = TRUE), top_n)
adt_top <- adt_avg[top_idx, ]

### Z-score by gene
adt_scaled <- t(scale(t(adt_top)))         # row-scale
adt_scaled <- adt_scaled[complete.cases(adt_scaled), ]  # drop rows with all NA

### Draw the heat-map
heat <- pheatmap(
  mat            = adt_scaled,
  color          = colorRampPalette(c("lightgrey", "blue", "red"))(100),
  fontsize_row   = 7,
  fontsize_col   = 10,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  main           = "Scaled ADT expression by annotated cluster",
  angle_col      = 45,
  silent         = TRUE
)

saveRDS(merged_seurat_log.harmony, file = "processed_seurat_object_MDS.rds")



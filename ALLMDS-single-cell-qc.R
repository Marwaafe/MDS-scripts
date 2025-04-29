An object of class Seurat 
36738 features across 531131 samples within 2 assays 
Active assay: RNA (36601 features, 0 variable features)
 1 layer present: counts
 1 other assay present: ADT
Error in `FetchData()`:
! None of the requested variables were found: 
Run `rlang::last_trace()` to see where the error occurred.
> rlang::last_trace()
<error/varsNotFoundError>
Error in `FetchData()`:
! None of the requested variables were found: 
---
Backtrace:
    ▆
 1. ├─base::subset(...)
 2. └─SeuratObject:::subset.Seurat(...)
 3.   ├─SeuratObject::WhichCells(...)
 4.   └─SeuratObject:::WhichCells.Seurat(...)
 5.     ├─SeuratObject::FetchData(...)
 6.     └─SeuratObject:::FetchData.Seurat(...)
Run rlang::last_trace(drop = FALSE) to see 1 hidden frame.
> rlang::last_trace(drop= FALSE)
<error/varsNotFoundError>
Error in `FetchData()`:
! None of the requested variables were found: 
---
Backtrace:
    ▆
 1. ├─base::subset(...)
 2. └─SeuratObject:::subset.Seurat(...)
 3.   ├─SeuratObject::WhichCells(...)
 4.   └─SeuratObject:::WhichCells.Seurat(...)
 5.     ├─SeuratObject::FetchData(...)
 6.     └─SeuratObject:::FetchData.Seurat(...)
 7.       └─rlang::abort(...)
> 

# Step 0: Load libraries
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

# Base folder path
base_path <- "/trinity/home/mafechkar/MDS_OUTS_CellRangerCount_9.0"

# Loop through each sample
for (sample_name in sample_names) {
  message("Processing sample: ", sample_name)
  
  # Set file path
  data_path <- file.path(
    base_path,
    paste0(sample_name, "_count_output"),
    paste0(sample_name, "_count"),
    "outs",
    "raw_feature_bc_matrix.h5"
  )
  
  # Step 2: Read data
  data_list <- Read10X_h5(data_path)
  
  # Step 3: Create Seurat object and add ADT
  seurat_obj <- CreateSeuratObject(counts = data_list[["Gene Expression"]])
  seurat_obj[["ADT"]] <- CreateAssayObject(counts = data_list[["Antibody Capture"]])
  
  # Step 4: Print summary
  print(seurat_obj)
  
  # Step 5: Calculate percent mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Step 6: Remove cells with missing QC values
  seurat_obj <- subset(seurat_obj, subset = !is.na(nFeature_RNA) & !is.na(nCount_RNA) & !is.na(percent.mt))
  
  # Step 7: Normalize RNA
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", assay = "RNA")
  
  # Step 8: Normalize ADT
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")
  
  # Step 9: Plot RNA QC violin plots (before filtering)
  p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
    geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2500, linetype = "dashed", color = "red")
  p2 <- VlnPlot(seurat_obj, features = "nCount_RNA")
  p3 <- VlnPlot(seurat_obj, features = "percent.mt")
  
  # Save violin plots before filtering
  qc_plot_before <- p1 + p2 + p3
  ggsave(filename = paste0(sample_name, "_QC_violin_before_filtering.png"), plot = qc_plot_before, width = 12, height = 5)

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
  
  # Step 12: Filter cells
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)
  
  # Step 13: Print number of cells after filtering
  cat("Cells after filtering:", ncol(seurat_obj), "\n")
  
  # Step 14: Replot QC violin plots after filtering
  p1_filtered <- VlnPlot(seurat_obj, features = "nFeature_RNA") +
    geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2500, linetype = "dashed", color = "red") +
    ggtitle("nFeature_RNA after filtering")
  p2_filtered <- VlnPlot(seurat_obj, features = "nCount_RNA") +
    ggtitle("nCount_RNA after filtering")
  p3_filtered <- VlnPlot(seurat_obj, features = "percent.mt") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    ggtitle("percent.mt after filtering")
  
  qc_plot_after <- p1_filtered + p2_filtered + p3_filtered
  
  # Save violin plots after filtering
  ggsave(filename = paste0(sample_name, "_QC_violin_after_filtering.png"), plot = qc_plot_after, width = 12, height = 5)

  # Step 15: Save filtered Seurat object
  saveRDS(seurat_obj, file = paste0(sample_name, "_filtered_normalized.rds"))
  
  message("Finished processing: ", sample_name, "\n")
}

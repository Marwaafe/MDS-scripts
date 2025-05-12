print(seurat_obj)
An object of class Seurat 
16915 features across 15911 samples within 2 assays 
Active assay: RNA (16778 features, 0 variable features)
 1 layer present: counts
 1 other assay present: ADT

////

> DefaultAssay(seurat_obj) <- 'RNA'
> rownames(seurat_obj)[1:10]
 [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"  "AL627309.2"  "AL627309.5"  "AL627309.4" 
 [9] "AP006222.2"  "AL732372.1" 
> DefaultAssay(seurat_obj) <- 'ADT'
> DefaultAssay(seurat_obj) <- 'ADT'
> rownames(seurat_obj)[1:10]
 [1] "CD86.1"     "CD274.1"    "TNFRSF14.1" "PVR.1"      "NECTIN2.1"  "CD47.1"     "CD48.1"     "CD40.1"     "CD40LG.1"  
[10] "CD52.1"    
> rownames(seurat_obj)
  [1] "CD86.1"      "CD274.1"     "TNFRSF14.1"  "PVR.1"       "NECTIN2.1"   "CD47.1"      "CD48.1"      "CD40.1"     
  [9] "CD40LG.1"    "CD52.1"      "CD3D.1"      "CD8A.1"      "NCAM1.1"     "CD19.1"      "CD33.1"      "ITGAX.1"    
 [17] "HLA-A.1"     "PTPRC.1"     "IL3RA.1"     "CD7.1"       "ENG.1"       "ITGA6.1"     "CCR4.1"      "CD4.1"      
 [25] "CD44.1"      "CD14.1"      "FCGR3A.1"    "IL2RA.1"     "PTPRC.2"     "PDCD1.1"     "TIGIT.1"     "nan"        
 [33] "nan.1"       "nan.2"       "nan.3"       "MS4A1.1"     "NCR1.1"      "PECAM1.1"    "MCAM.1"      "IGHM.1"     
 [41] "CD5.1"       "CXCR3.1"     "CCR5.1"      "FCGR2A.1"    "CCR6.1"      "CXCR5.1"     "ITGAE.1"     "CD69.1"     
 [49] "SELL.1"      "KLRB1.1"     "CTLA4.1"     "LAG3.1"      "KLRG1.1"     "CD27.1"      "LAMP1.1"     "FAS.1"      
 [57] "TNFRSF4.1"   "HLA-DRA.1"   "CD1C.1"      "ITGAM.1"     "FCGR1A.1"    "THBD.1"      "CD1D.1"      "KLRK1.1"    
 [65] "CR1.1"       "B3GAT1.1"    "BTLA.1"      "ICOS.1"      "CD58.1"      "ENTPD1.1"    "CX3CR1.1"    "CD24.1"     
 [73] "CR2.1"       "ITGAL.1"     "CD79B.1"     "CD244.1"     "SIGLEC1.1"   "ITGB7.1"     "TNFRSF13C.1" "GP1BB.1"    
 [81] "ICAM1.1"     "SELP.1"      "IFNGR1.1"    "nan.4"       "nan.5"       "nan.6"       "nan.7"       "IL2RB.1"    
 [89] "TNFRSF13B.1" "FCER1A.1"    "ITGA2B.1"    "TNFRSF9.1"   "CD163.1"     "CD83.1"      "IL4R.1"      "ANPEP.1"    
 [97] "CD2.1"       "CD226.1"     "ITGB1.1"     "CLEC4C.1"    "ITGA2.1"     "CD81.1"      "IGHD.1"      "ITGB2.1"    
[105] "CD28.1"      "CD38.1"      "IL7R.1"      "PTPRC.3"     "CD22.1"      "TFRC.1"      "DPP4.1"      "CD36.1"     
[113] "KIR2DL1.1"   "ITGA1.1"     "ITGA4.1"     "NT5E.1"      "nan.8"       "nan.9"       "OLR1.1"      "KIR2DL3.1"  
[121] "KIR3DL1.1"   "SLAMF7.1"    "CD99.1"      "CLEC12A.1"   "SLAMF6.1"    "KLRD1.1"     "IGKC.1"      "LILRB1.1"   
[129] "FCER2.1"     "nan.10"      "SIGLEC7.1"   "ADGRG1.1"    "HLA-E.1"     "CD82.1"      "CD101.1"     "C5AR1.1"    
[137] "GGT1.2"  


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

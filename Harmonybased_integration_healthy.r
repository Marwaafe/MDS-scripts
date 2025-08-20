library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)

sample_names_healthy <- c(
  "HealthyBM1", "HealthyBM2" , "HealthyBM3" , "HealthyBM4" , "HealthyBM5" , "HealthyBM6", "HealthyBM7"
  , "HealthyBM8" , "HealthyBM9"
)

seurat_list_Healthy <- lapply(sample_names_healthy, function(name) {
  readRDS(file.path("/trinity/home/mafechkar/Healthy_BoneM", paste0(name, "_SeuratObj.rds")))
})
names(seurat_list_Healthy) <- sample_names_healthy
merged_Healthy <- merge(
  x = seurat_list_Healthy[[1]],
  y = seurat_list_Healthy[-1],
  add.cell.ids = sample_names_healthy,
  project = "Healthy_Merged",
  merge.data = FALSE  # keep normalized RNA and ADT layers
)

merged_Healthy$orig.ident <- sapply(strsplit(colnames(merged_Healthy), "_"), `[`, 1)


merged_Healthy
head(VariableFeatures(merged_Healthy))

merged_Healthy <- subset(merged_Healthy, 
                         subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 5000 & 
                           percent.mt < 10)


#Normalize using LogNormalize with median scaling factor
median_umi <- median(merged_Healthy$nCount_RNA)
print(median_umi)
merged_Healthy_log <- NormalizeData(merged_Healthy, normalization.method = "LogNormalize", scale.factor = median_umi)
str(merged_Healthy_log)
merged_Healthy_log@meta.data
head(merged_Healthy_log[[]])
VlnPlot(merged_Healthy_log, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



#Detecting Batch effect
merged_Healthy_log <- FindVariableFeatures(merged_Healthy_log)
merged_Healthy_log <- ScaleData(merged_Healthy_log)
merged_Healthy_log <- RunPCA(merged_Healthy_log)
ElbowPlot(merged_Healthy_log)
merged_Healthy_log <- RunUMAP(merged_Healthy_log, dim = 1:20, reduction = 'pca')
distinct_colors <- c(
  "#E31A1C",  # red
  "#1F78B4",  # strong blue
  "#33A02C",  # green
  "#FF7F00",  # vivid orange
  "#6A3D9A",  # purple
  "#A6CEE3",  # light blue
  "#B2DF8A",  # light green
  "#FB9A99",  # light red/pink
  "#CAB2D6",  # lavender
  "#FDBF6F",  # peach
  "#FFED6F",  # yellow
  "#B15928",  # brown
  "#66C2A5",  # teal
  "#882255",  # salmon
  "#F0027F",  # steel blue
  "#E78AC3",  # pale pink
  "#117733",  # lime
  "#4169E1",  # warm yellow 
  "#00CED1",  # dark turquoise
  "#DC143C",  # crimson
  "#8B4513",  # saddle brown
  "#FFD700",  # gold
  "#7FFF00"   # chartreuse green
)

before_Healthy<- DimPlot(merged_Healthy_log, reduction = "umap", group.by = "orig.ident") +
  scale_color_manual(values = distinct_colors) +
  ggtitle("Healthy Donors Colored by Samples") +
  theme(aspect.ratio = 1)
# run Harmony -----------
merged_Healthy_log.harmony <- merged_Healthy_log %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = FALSE)

merged_Healthy_log.harmony@reductions

merged_Healthy_log.harmony.embed <- Embeddings(merged_Healthy_log.harmony, "harmony")
merged_Healthy_log.harmony.embed[1:10,1:10]


# UMAP and clustering using ** Harmony embeddings instead of PCA **
merged_Healthy_log.harmony <- merged_Healthy_log.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after_Healthy<- DimPlot(merged_Healthy_log.harmony, reduction = 'umap', group.by = 'orig.ident') +
  scale_color_manual(values = distinct_colors) +
  ggtitle("After HARMONY batch correction") +
  theme(aspect.ratio = 1)

before_Healthy|after_Healthy

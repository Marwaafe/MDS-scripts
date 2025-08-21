

set.seed(1234)

library(Seurat)
install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("metap")
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

install.packages('harmony')
install.packages('Rcpp')
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)

sample_names <- c(
  "MDS001-09-203", "MDS005-09-247", "MDS006-08-249", "MDS010-09-299",
  "MDS016-09-478", "MDS023-10-053", "MDS029-10-118", "MDS038-10-241",
  "MDS059-10-531", "MDS065-10-609", "MDS154-13-486", "MDS155-13-606",
  "MDS167-13-913", "MDS169-13-919", "MDS180-14-164", "MDS189-14-527",
  "MDS201-15-093", "MDS212-15-463"
)

seurat_list <- lapply(sample_names, function(name) {
  readRDS(file.path("/trinity/home/mafechkar/MDS_Data", paste0(name, "_SeuratObj.rds")))
})
names(seurat_list) <- sample_names
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_names,
  project = "MDS_Merged",
  merge.data = FALSE  # keep normalized RNA and ADT layers
)

merged_seurat$orig.ident <- sapply(strsplit(colnames(merged_seurat), "_"), `[`, 1)

#Normalize using LogNormalize with median scaling factor
median_umi <- median(merged_seurat$nCount_RNA)
print(median_umi)
merged_seurat_log <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = median_umi)
str(merged_seurat_log)
merged_seurat_log@meta.data

#Detecting Batch effect
merged_seurat_log <- FindVariableFeatures(merged_seurat_log)
merged_seurat_log <- ScaleData(merged_seurat_log)
merged_seurat_log <- RunPCA(merged_seurat_log)
ElbowPlot(merged_seurat_log)
merged_seurat_log <- RunUMAP(merged_seurat_log, dim = 1:20, reduction = 'pca')
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
  "#4169E1",  # royal blue
  "#00CED1",  # Dark Turquoise
  "#FFD700",  # Gold
  "#7FFF00",  # Chartreuse
  "#DC143C",  # Crimson
  "#FF69B4",  # Hot Pink
  "#20B2AA",  # Light Sea Green
  "#8B008B",  # Dark Magenta
  "#00FA9A",  # Medium Spring Green
  "#FF4500",  # Orange Red
  "#DAA520",  # Goldenrod
  "#5F9EA0",  # Cadet Blue
  "#BA55D3",  # Medium Orchid
  "#ADFF2F",  # Green Yellow
  "#2E8B57",  # Sea Green
  "#FF1493",  # Deep Pink
  "#8B4513",  # Saddle Brown
  "#40E0D0",  # Turquoise
  "#9932CC",  # Dark Orchid
  "#A52A2A",  # Brick Red
  "#DC75CD",  # Orchid Pink
  "#228B22",  # Forest Green
  "#FF6347",  # Tomato Red
  "#1E90FF",  # Dodger Blue
  "#FFDAB9",  # Peach Puff
  "#000000",  # Black
  "#808000",  # Olive
  "#FF00FF",  # Fuchsia
  "#008080",  # Teal
  "#C71585",  # Medium Violet Red
  "#191970",  # Midnight Blue
  "#00BFFF",  # Deep Sky Blue
  "#006400",  # Dark Green
  "#FF8C00",  # Dark Orange
  "#483D8B",  # Dark Slate Blue
  "#B22222",  # Firebrick
  "#9ACD32",  # Yellow Green
  "#8FBC8F"   # Dark Sea Green
)


before_MDS <-DimPlot(merged_seurat_log, reduction = "umap", group.by = "orig.ident") +
  scale_color_manual(values = distinct_colors) +
  ggtitle("MDS colored by samples") +
  theme(aspect.ratio = 1)

DimPlot(merged_seurat_log, reduction = "umap", group.by = "orig.ident")

# run Harmony -----------
merged_seurat_log.harmony <- merged_seurat_log %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = FALSE)

merged_seurat_log.harmony@reductions

merged_seurat_log.harmony.embed <- Embeddings(merged_seurat_log.harmony, "harmony")
merged_seurat_log.harmony.embed[1:10,1:10]


# UMAP and clustering using ** Harmony embeddings instead of PCA **
merged_seurat_log.harmony <- merged_seurat_log.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8)

# visualize 
after_MDS <-DimPlot(merged_seurat_log.harmony, reduction = 'umap', group.by = 'orig.ident') +
  scale_color_manual(values = distinct_colors) +
  ggtitle("After HARMONY batch correction") +
  theme(aspect.ratio = 1)

before_MDS|after_MDS

# Run this BEFORE FindAllMarkers
merged_seurat_log.harmony <- JoinLayers(merged_seurat_log.harmony)
install.packages('devtools')
devtools::install_github('immunogenomics/presto', force = TRUE)
library(presto)

markers <- FindAllMarkers(
  merged_seurat_log.harmony,
  min.pct = 0.2,
  min.diff.pct = 0.2,
  only.pos = TRUE,
  verbose = TRUE,
  slot = "data"
) 

# Then run marker detection
markers <- FindAllMarkers(
  merged_seurat_log.harmony,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)


top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

top_markers
print(top_markers, n = 115)

# Top 5 genes per cluster
top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))


# Ensure uniqueness of gene names (DotPlot can't handle duplicates in a named list)
top5_genes_ordered <- unique(top5_markers$gene)

# Plot using DotPlot
DotPlot(
  merged_seurat_log.harmony,
  features = top5_genes_ordered,
  group.by = "seurat_clusters"
) + RotatedAxis()

VlnPlot(
  merged_seurat_log.harmony,
  features = top5_markers$gene,
  group.by = "seurat_clusters",
  stack = TRUE,
  flip = TRUE,
  pt.size = 0
) + NoLegend()

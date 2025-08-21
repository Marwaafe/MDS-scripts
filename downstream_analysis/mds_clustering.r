# script to identify cluster identity -----------------
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions

# setwd("~/")

set.seed(1234)

library(Seurat)
install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("metap")
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

# load data
merged_seurat_log.harmony <- readRDS(merged_seurat_log.harmony)
str(merged_seurat_log.harmony)
View(merged_seurat_log.harmony@meta.data)

library(clustree)

merged_seurat_log.harmony <- FindNeighbors(merged_seurat_log.harmony, reduction = 'harmony', dims = 1:20)  #  can tweak dims

merged_seurat_log.harmony <- FindClusters(merged_seurat_log.harmony, resolution = seq(from = 0.1, to = 0.9, by = 0.1))
merged_seurat_log.harmony@meta.data
clust_seurat <- merged_seurat_log.harmony@meta.data %>% dplyr::select(dplyr::contains("RNA_snn_res."))
clustree(clust_seurat, prefix="RNA_snn_res.")


# 5. re-run clustering at that specific resolution
merged_seurat_log.harmony <- FindClusters(merged_seurat_log.harmony, resolution = 0.8)  # resolution controls granularity

# 6. Run UMAP
merged_seurat_log.harmony <- RunUMAP(merged_seurat_log.harmony, reduction = 'harmony', dims = 1:20)


# visualize data
clusters <- DimPlot(merged_seurat_log.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE) + scale_color_manual(values = distinct_colors) +
  ggtitle("MDS colored by clusers 0.8") +
  theme(aspect.ratio = 1)
condition <- DimPlot(merged_seurat_log.harmony, reduction = "umap", group.by = "orig.ident") +
  scale_color_manual(values = distinct_colors) +
  ggtitle("MDS colored by samples") +
  theme(aspect.ratio = 1)

condition|clusters

# findAll markers -----------------

FindAllMarkers(merged_seurat_log.harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

markers <- FindAllMarkers(
  merged_seurat_log.harmony,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",   # default
  slot = "data"          # default slot for normalized data
)

# BEFORE FindAllMarkers
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



#Cluster 17 quality count check

# Create a vector of gene counts for all cells (NA for non-cluster 17)
gene_counts_17 <- rep(NA, ncol(merged_seurat_log.harmony))
names(gene_counts_17) <- colnames(merged_seurat_log.harmony)
cluster_17_cells <- WhichCells(merged_seurat_log.harmony, idents = "17")
gene_counts_17[cluster_17_cells] <- merged_seurat_log.harmony$nFeature_RNA[cluster_17_cells]

# Add this vector as a new meta.data column
merged_seurat_log.harmony$gene_count_cluster17 <- gene_counts_17

# Plot using FeaturePlot (no na.col)
FeaturePlot(
  merged_seurat_log.harmony,
  features = "gene_count_cluster17",
  cols = c("lightgrey", "gold", "red"),
  order = TRUE
)
FeaturePlot(
  merged_seurat_log.harmony,
  features = "gene_count_cluster17",
  cols = viridis::viridis(100),
  order = TRUE
)
# Create vector of NA values for all cells
nCount_cluster17 <- rep(NA, ncol(merged_seurat_log.harmony))
names(nCount_cluster17) <- colnames(merged_seurat_log.harmony)

# Fill in values for cluster 17 cells
cluster17_cells <- WhichCells(merged_seurat_log.harmony, idents = "17")
nCount_cluster17[cluster17_cells] <- merged_seurat_log.harmony$nCount_RNA[cluster17_cells]

# Add to metadata
merged_seurat_log.harmony$nCount_cluster17 <- nCount_cluster17

# Plot
FeaturePlot(
  merged_seurat_log.harmony,
  features = "nCount_cluster17",
  cols = c("lightgrey", "yellow", "red"),
  pt.size = 0.5,
  order = TRUE
)



# Subset cluster 17 from the full object
cells_18 <- WhichCells(merged_seurat_log.harmony, idents = 18)

# Create a named vector of nFeature_RNA values just for cluster 17 cells
gene_counts <- merged_seurat_log.harmony@meta.data[cells_18, "nFeature_RNA"]

# Add a metadata column with NA for all, then fill only for cluster 17
merged_seurat_log.harmony$cluster18_gene_counts <- NA
merged_seurat_log.harmony$cluster18_gene_counts[cells_18] <- gene_counts

# Plot UMAP and color only cluster 17 by gene counts
FeaturePlot(
  merged_seurat_log.harmony,
  features = "cluster18_gene_counts",
  reduction = "umap",
  cols = c("yellow", "blue"),
  pt.size = 0.3
) + 
  ggtitle("UMAP: genes per cell for cluster 18")


#Cluster 18 quality count check
cluster18 <- subset(merged_seurat_log.harmony, idents = c(18))

gene_counts <- cluster18@meta.data$nFeature_RNA

plot(
  sort(gene_counts),
  main = "genes per cell in cluster 18",
  xlab = "cell",
  ylab = "genes_per_cell",
  pch = 1,
  col = "black"
)

cluster18[["percent.mt"]] <- PercentageFeatureSet(cluster18, pattern = "^MT-")

# Plot ordered percent.mt
mt_percent <- cluster18@meta.data$percent.mt
plot(
  sort(mt_percent),
  main = "percent.mt in cluster 18",
  xlab = "cell",
  ylab = "% mitochondrial genes",
  pch = 1,
  col = "black"
)

# findConserved markers -------------

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(merged_seurat_log) <- 'RNA'

DefaultAssay(merged_seurat_log.harmony)

markers_cluster0 <- FindConservedMarkers(merged_seurat_log.harmony,
                                         ident.1 = 0,
                                         grouping.var = 'orig.ident')

head(markers_cluster0)
str(markers_cluster0)

# let's visualize top features
FeaturePlot(merged_seurat_log.harmony, features = c('LTB'), min.cutoff = 'q10')


# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))





# rename cluster 3 ident
Idents(merged_seurat_log.harmony)
merged_seurat_log.harmony <- RenameIdents(merged_seurat_log.harmony, `3` = 'CD16 Mono')

DimPlot(merged_seurat_log.harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
View(merged_seurat_log.harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(merged_seurat_log.harmony) <- merged_seurat_log.harmony@meta.data$seurat_annotations
Idents(merged_seurat_log.harmony)

DimPlot(merged_seurat_log.harmony, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
merged_seurat_log.harmony$celltype.cnd <- paste0(merged_seurat_log.harmony$seurat_annotations,'_', merged_seurat_log$stim)
View(merged_seurat_log.harmony@meta.data)
Idents(merged_seurat_log.harmony) <- merged_seurat_log.harmony$celltype.cnd

DimPlot(merged_seurat_log.harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(merged_seurat_log.harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(merged_seurat_log.harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')

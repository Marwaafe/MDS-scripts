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
merged_Healthy_log.harmony <- readRDS(merged_Healthy_log.harmony)
str(merged_Healthy_log.harmony)
View(merged_Healthy_log.harmony@meta.data)

library(clustree)

merged_Healthy_log.harmony <- FindNeighbors(merged_Healthy_log.harmony, reduction = 'harmony', dims = 1:20)  #  can tweak dims

merged_Healthy_log.harmony <- FindClusters(merged_Healthy_log.harmony, resolution = seq(from = 0.1, to = 0.9, by = 0.1))
merged_Healthy_log.harmony@meta.data
clust_healthy <- merged_Healthy_log.harmony@meta.data %>% dplyr::select(dplyr::contains("RNA_snn_res."))
clustree(clust_healthy, prefix="RNA_snn_res.")


# 5. re-run clustering at that specific resolution
merged_Healthy_log.harmony <- FindClusters(merged_Healthy_log.harmony, resolution = 0.8)  # resolution controls granularity

# 6. Run UMAP
merged_Healthy_log.harmony <- RunUMAP(merged_Healthy_log.harmony, reduction = 'harmony', dims = 1:20)


# visualize data
clusters_healthy <- DimPlot(merged_Healthy_log.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE) + scale_color_manual(values = distinct_colors) +
  ggtitle("Healthy colored by clusters 0.8") +
  theme(aspect.ratio = 1)
condition_healthy <- DimPlot(merged_Healthy_log.harmony, reduction = "umap", group.by = "orig.ident") +
  scale_color_manual(values = distinct_colors) +
  ggtitle("Healthy colored by samples") +
  theme(aspect.ratio = 1)

condition_healthy|clusters_healthy

# findAll markers -----------------

FindAllMarkers(merged_Healthy_log.harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

markers_healthy <- FindAllMarkers(
  merged_Healthy_log.harmony,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",   # default
  slot = "data"          # default slot for normalized data
)

# Run this BEFORE FindAllMarkers
merged_Healthy_log.harmony <- JoinLayers(merged_Healthy_log.harmony)
install.packages('devtools')
devtools::install_github('immunogenomics/presto', force = TRUE)
library(presto)

# Then run marker detection

markers_healthy <- FindAllMarkers(
  merged_Healthy_log.harmony,
  min.pct = 0.2,
  min.diff.pct = 0.2,
  only.pos = TRUE,
  verbose = TRUE,
  slot = "data"
) 

top_markers_healthy <- markers_healthy %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

top_markers_healthy
print(top_markers_healthy, n = 23)

FeaturePlot(merged_Healthy_log.harmony, features = c("CD3E", "CD4", "CD8A"))
FeaturePlot(merged_Healthy_log.harmony, features = c("MS4A1", "CD79A"))
FeaturePlot(merged_Healthy_log.harmony, features = c("LYZ", "CD14", "FCGR3A"))
FeaturePlot(merged_Healthy_log.harmony, features = c("NKG7", "GNLY"))
FeaturePlot(merged_Healthy_log.harmony, features = c("MPO", "PRTN3"))
FeaturePlot(merged_Healthy_log.harmony, features = c("SDC1", "IGLC2"))


DimPlot(merged_Healthy_log.harmony, group.by = "seurat_clusters", label = TRUE)

# findConserved markers -------------

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(merged_seurat_log) <- 'RNA'

DefaultAssay(merged_Healthy_log.harmony)

markers_Healthy_cluster0 <- FindConservedMarkers(merged_Healthy_log.harmony,
                                         ident.1 = 0,
                                         grouping.var = 'orig.ident')

head(markers_Healthy_cluster0)
str(markers_Healthy_cluster0)

# let's visualize top features
FeaturePlot(merged_Healthy_log.harmony, features = c('GNLY'), min.cutoff = 'q10')


# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))





# rename cluster 3 ident
Idents(merged_Healthy_log.harmony)
merged_Healthy_log.harmony <- RenameIdents(merged_Healthy_log.harmony, `4` = 'Memory T-Cell')

DimPlot(merged_Healthy_log.harmony, reduction = 'umap',  label = T)

# cells already have annotations provided in the metadata
View(merged_Healthy_log.harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(merged_Healthy_log.harmony) <- merged_Healthy_log.harmony@meta.data$seurat_annotations
Idents(merged_healthy_log.harmony)

DimPlot(merged_healthy_log.harmony, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
merged_Healthy_log.harmony$celltype.cnd <- paste0(merged_Healthy_log.harmony$seurat_annotations,'_', merged_seurat_log$stim)
View(merged_Healthy_log.harmony@meta.data)
Idents(merged_Healthy_log.harmony) <- merged_Healthy_log.harmony$celltype.cnd

DimPlot(merged_Healthy_log.harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(merged_Healthy_log.harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_Healthy_cluster3)


FeaturePlot(merged_Healthy_log.harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')

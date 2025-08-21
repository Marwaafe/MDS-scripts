sample_names <- c(
  "HealthyBM1", "HealthyBM2" , "HealthyBM3" , "HealthyBM4" , "HealthyBM5" , "HealthyBM6", "HealthyBM7"
  , "HealthyBM8" , "HealthyBM9"
)

seurat_list <- lapply(sample_names, function(name) {
  readRDS(file.path("/trinity/home/mafechkar/Healthy_BoneM", paste0(name, "_SeuratObj.rds")))
})
names(seurat_list) <- sample_names
merged_seu_HD <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_names,
  project = "Healthy_Merged",
  merge.data = TRUE  # keep normalized RNA and ADT layers
)

merged_seu_HD$orig.ident <- sapply(strsplit(colnames(merged_seu_HD), "_"), `[`, 1)


merged_seu_HD
head(VariableFeatures(merged_seu_HD))

merged_seu_HD <- NormalizeData(merged_seu_HD)
merged_seu_HD <- FindVariableFeatures(merged_seu_HD)
merged_seu_HD <- ScaleData(merged_seu_HD)
merged_seu_HD <- RunPCA(merged_seu_HD)
merged_seu_HD <- FindNeighbors(merged_seu_HD, dims = 1:30, reduction = "pca")
merged_seu_HD <- FindClusters(merged_seu_HD, resolution = 0.5, cluster.name = "unintegrated_clusters")
merged_seu_HD <- RunUMAP(merged_seu_HD, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(merged_seu_HD, reduction = "umap.unintegrated", group.by = "orig.ident")

DimPlot(merged_seu_HD, reduction = "umap.unintegrated", group.by = "seurat_clusters", label = FALSE)




saveRDS(merged_seu_HD, file = "/trinity/home/mafechkar/MDS_Data/Healthy_Merged_SeuratObj.rds")

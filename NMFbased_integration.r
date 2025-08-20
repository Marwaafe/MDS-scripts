library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(UCell)
library(RColorBrewer)
library(GeneNMF)

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

#Detecting Batch effect
merged_seurat_log <- FindVariableFeatures(merged_seurat_log)
merged_seurat_log <- ScaleData(merged_seurat_log)
merged_seurat_log <- RunPCA(merged_seurat_log)
ElbowPlot(merged_seurat_log)
merged_seurat_log <- RunUMAP(merged_seurat_log, dim = 1:20, reduction = 'pca')
before_NMF <-DimPlot(merged_seurat_log, reduction = "umap", group.by = "orig.ident")

#set RNA assay as default and split object by sample 
DefaultAssay(merged_seurat_log) <- "RNA"
seurat_list <- SplitObject(merged_seurat_log, split.by = "orig.ident")

common_genes <- Reduce(intersect, lapply(seurat_list, function(x) rownames(x)))

seurat_list <- lapply(seurat_list, function(x) {
  x[common_genes,]
})

#Run NMF across multiple k values
geneNMF.programs <- multiNMF(seurat_list, assay = "RNA", k = 12:20, min.exp = 0.01)



#Collapse those programs into consensus metaprograms (MPs)
geneNMF.metaprograms <- getMetaPrograms(
  geneNMF.programs,
  metric = "cosine",
  weight.explained = 0.5,
  nMP = 18
)

ph <- plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0.1, 1))

geneNMF.metaprograms$metaprograms.metrics

keep_mps <- c("MP1", "MP2", "MP4", "MP5", "MP6", "MP7", "MP9", "MP8", "MP10", "MP11", "MP15", "MP16", "MP17")
filtered_metaprograms <- geneNMF.metaprograms$metaprograms.genes[keep_mps]


str(filtered_metaprograms)

#identify non wanted genes by prefix 
ALs <- grep('^AL[0-9]',rownames(merged_seurat_log[['RNA']]),value = TRUE)
ACs <- grep('^AC[0-9]',rownames(merged_seurat_log[['RNA']]),value = TRUE)
ADs <- grep('^AD[0-9]',rownames(merged_seurat_log[['RNA']]),value = TRUE)
APs <- grep('^AP[0-9]',rownames(merged_seurat_log[['RNA']]),value = TRUE)
EFs <-  grep('^EEF',rownames(merged_seurat_log[['RNA']]),value = TRUE)
LINCs <-  grep('^LINC',rownames(merged_seurat_log[['RNA']]),value = TRUE)

#combine into a single gene list
genes_to_remove <- c('MALAT1','LNCAROD','IGKC', ALs, ACs, ADs,APs,EFs,LINCs)
genes_to_keep <- setdiff(rownames(merged_seurat_log), genes_to_remove)
merged_seurat_log <- subset(merged_seurat_log, features = genes_to_keep)

#rerun metaprogram analysis on clean data
geneNMF.metaprograms <- getMetaPrograms(
  geneNMF.programs,
  metric = "cosine",
  weight.explained = 0.5,
  nMP = 18
)

ph <- plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0.1, 1))
geneNMF.metaprograms$metaprograms.metrics
keep_mps <- c("MP1", "MP2", "MP4", "MP5", "MP6", "MP7", "MP9", "MP8", "MP10", "MP11", "MP15", "MP16", "MP17")
filtered_metaprograms <- geneNMF.metaprograms$metaprograms.genes[keep_mps]

#interpretation of each gene prgram 
lapply(geneNMF.metaprograms$metaprograms.genes, head)
lapply(geneNMF.metaprograms$metaprograms.genes.weights, head)

geneNMF.metaprograms$metaprograms.genes.weights$MP1

#Functional Annotation via GSEA
install.packages("msigdbr")
install.packages("BiocManager")
BiocManager::install("fgsea")
library(msigdbr)
library(fgsea)

top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(merged_seurat_log), category = "C5", subcategory = "GO:BP")
})
head(top_p$MP1) #Loop through all to interpret the program sepearetly 
head(top_p$MP2)
head(top_p$MP4)
head(top_p$MP5)
head(top_p$MP6)
head(top_p$MP7)
head(top_p$MP8)
head(top_p$MP9)
head(top_p$MP10)
head(top_p$MP11)
head(top_p$MP15)
head(top_p$MP16)
head(top_p$MP17)
#score each cell using UCell
#keeping only th metaprograms I'm using
mp.genes <- geneNMF.metaprograms$metaprograms.genes[c("MP1", "MP2", "MP4", "MP5", "MP6", "MP7", "MP9", "MP8", "MP10", "MP11", "MP16")]

merged_seurat_log_NMF<- merged_seurat_log

#score each cell
merged_seurat_log_NMF <- UCell::AddModuleScore_UCell(
  merged_seurat_log_NMF,
  features = mp.genes,
  ncores = 4,
  name = ""
)

#Visualizing metaprogram activity across samples
VlnPlot(merged_seurat_log_NMF, features = names(mp.genes), group.by = "orig.ident", pt.size = 0, ncol = 4)

#Extract MP scores from metadata
#matrix of Mp scores already stored in metadata by UCell
matrix <- merged_seurat_log_NMF@meta.data[, names(mp.genes)]
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_", seq_len(ncol(dimred)))

#create new dimreduc slot
merged_seurat_log_NMF@reductions[["MPsignatures"]] <- new(
  "DimReduc",
  cell.embeddings = dimred,
  assay.used = "RNA",
  key = "MP_",
  global = FALSE
)
merged_seurat_log@reductions$NMF


#run umap on MP signature space 
set.seed(123)
merged_seurat_log_NMF <- RunUMAP(
  merged_seurat_log_NMF,
  reduction = "MPsignatures",
  dims = 1:ncol(dimred),
  metric = "euclidean",
  reduction.name = "umap_MP"
)

#visualize the new integrated umap
after_NMF <- DimPlot(merged_seurat_log_NMF, reduction = "umap_MP", group.by = "orig.ident") +
  theme(aspect.ratio = 1)

before_NMF | after_NMF

# Generate UMAP with distinct colors
DimPlot(merged_seurat_log_NMF, reduction = "umap_MP", group.by = "orig.ident") +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#14168e",
    "#8DA0CB", "#8e1440", "#A6D854", "#FFD92F", "#E5C494",
    "#B3B3B3", "#1B9E77", "#0ad9e3", "#7570B3", "#E7298A"
  )) +
  ggtitle("UMAP colored by sample (orig.ident)") +
  theme(aspect.ratio = 1)


#show MP activity across UMAP
library(viridis)
FeaturePlot(merged_seurat_log_NMF, features = names(mp.genes), reduction = "umap_MP", ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

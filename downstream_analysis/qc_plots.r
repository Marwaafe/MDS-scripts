#QC
counts <- seu_mds1@assays$RNA$counts
counts[1:10,1:3]
genes_per_cell <- Matrix::colSums(counts>0) # count a gene only if it has non-zero reads mapped.
counts_per_cell <- Matrix::colSums(counts)
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
MIN_GENES_PER_CELL <- 600
MAX_GENES_PER_CELL <- 6000
# now replot with the thresholds being shown:
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
abline(h=MIN_GENES_PER_CELL, col='magenta')  # lower threshold
abline(h=MAX_GENES_PER_CELL, col='gold') # upper threshold

#MT
seu_mds1[["percent.mt"]] <- PercentageFeatureSet(seu_mds1, pattern = "^MT-")
mito_genes <- grep("^mt-", rownames(counts) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(counts[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab =
       "percentage mitochondrial counts")

MAX_PCT_MITO <- 10
plot(sort(pct_mito))
abline(h=MAX_PCT_MITO, col='red')
seu_mds1 <- subset(seu_mds1, subset = nFeature_RNA >100 & nFeature_RNA <2000 & percent.mt<10)
dim(seu_mds1)
seu_mds1$IDs <- seu_mds1$orig.ident
levels(seu_mds1$IDs) <- 'HL_control'

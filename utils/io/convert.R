require(Seurat)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
name <- sub('.h5ad', '', filename, fixed = TRUE)
counts_tmp = paste0(name,"_counts.csv")
metadata_tmp = paste0(name,"_metadata.csv")


seu <- CreateSeuratObject(counts=t(read.csv(counts_tmp, row.names=1)), assay="RNA")
seu@meta.data <- read.csv(metadata_tmp, row.names=1)

saveRDS(seu, paste0(name,".rds"))

file.remove(counts_tmp)
file.remove(metadata_tmp)
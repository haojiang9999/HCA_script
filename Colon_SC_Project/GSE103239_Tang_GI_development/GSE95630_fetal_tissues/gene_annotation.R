library(mygene)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene", version = "3.8")
unknown <- geneNames[!geneNames %in% gencode.v19$gene_name]
saveRDS(unknown, file = "unknown genes.rds")
getGenes(unknown)

#### Samples annotation
colnames(Combine_df) <- gsub("_scTrioSeq2Rna_scTrioSeq2Rna_", "_scTrioSeq2Rna_", colnames(Combine_df))
Pheno.688 <- list()
Pheno.688$cellType <- sapply(strsplit(as.character(colnames(Combine_df)), "_"), "[[", 4 )
Pheno.688$cellType2 <-substr(Pheno.688$cellType, start = 1, stop = 2)
Pheno.688$pateint <- sapply(strsplit(as.character(colnames(Combine_df)), "_"), "[[", 3 )
Pheno.688$sampleNames <- colnames(Combine_df)
#### read lineage info
### Sad there were no cells matched cellSublineage info
cellSublineage <- read.csv("./Figure1_cell_sublineage_190710.csv")
sampleNamesShort <- sapply(strsplit(as.character(Pheno.688$sampleNames), "na_"), "[[", 2 )
sampleNamesShort <- sapply(strsplit(as.character(sampleNamesShort), ".TPM"), "[[", 1 )
sampleNamesShort %in% as.character(cellSublineage$Sample)
as.character(cellSublineage$Sample) %in% sampleNamesShort

saveRDS(Pheno.688, file = "Pheno.list.688.rds")

#### 3.convert ensembleID to symbol
#### 2.Test.R
COAD_UCSC_Toil_tpm_dataset <- readRDS("COAD_UCSC_Toil_tpm_dataset.rds")
COAD.RSEM.gene.tpm <- COAD_UCSC_Toil_tpm_dataset$COAD.RSEM.gene.tpm
COAD.RSEM.gene.tpm[1:5,1:5]
summary(colSums(COAD.RSEM.gene.tpm))
#head(gencode.v23.annotation)
## convert ensembleID to symbol 
geneMatch <- match(rownames(COAD_RSEM_gene_tpm),gencode.v23.annotation$V1)
geneSymbol <- as.character(gencode.v23.annotation[geneMatch,]$V2)
COAD_tpm_symbol <- COAD_RSEM_gene_tpm
rownames(COAD_tpm_symbol) <- geneSymbol
COAD_tpm_symbol[1:5,1:5]
COAD_UCSC_Toil_tpm_dataset$COAD.RSEM.gene.tpm.symbol <- COAD_tpm_symbol
saveRDS(COAD_UCSC_Toil_tpm_dataset, file = "COAD_UCSC_Toil_tpm_dataset.rds")

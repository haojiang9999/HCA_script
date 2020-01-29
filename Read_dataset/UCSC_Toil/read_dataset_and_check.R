#### read dataset and check
BRCA_UCSC_Toil_tpm_dataset <- readRDS("./BRCA_UCSC_Toil_tpm_dataset.rds")
BRCA.RSEM.gene.tpm <- BRCA_UCSC_Toil_tpm_dataset$BRCA.RSEM.gene.tpm
BRCA.RSEM.gene.tpm[1:5,1:5]
TCGA.sampleID <- as.character(phenotype_sp$sample)
TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
rownames(phenotype_sp) <- TCGA.sampleID
phenotype_sp["TCGA.BH.A0BQ.11",]

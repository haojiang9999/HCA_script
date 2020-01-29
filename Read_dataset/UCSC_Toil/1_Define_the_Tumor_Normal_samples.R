COAD_UCSC_Toil_tpm_dataset <- readRDS("COAD_UCSC_Toil_tpm_dataset.rds")
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset$COAD.pheno
### Seperate the Normal and Tumor samples
sampleID<- as.character(COAD.pheno$sample)
table(unlist(lapply(strsplit(sampleID,"-"), '[[', 4)))
sampleTypes <- unlist(lapply(strsplit(sampleID,"-"), '[[', 4))
sampleTypes[sampleTypes=="01"]<- "Tumor"
sampleTypes[sampleTypes=="11"]<- "Nomal"
COAD.pheno <- cbind(COAD.pheno , sampleTypes)

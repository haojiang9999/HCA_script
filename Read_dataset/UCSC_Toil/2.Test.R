#### 2.Test.R
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706 <- readRDS("COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
COAD.RSEM.gene.tpm <- COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.RSEM.gene.tpm
COAD.RSEM.gene.tpm[1:5,1:5]
summary(colSums(COAD.RSEM.gene.tpm))

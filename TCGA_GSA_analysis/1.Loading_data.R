#### 1.Loading_data.R
#1.read TCGA_COAD TPM data ####
COAD_UCSC_Toil_tpm_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset.rds")
COAD.RSEM.gene.tpm.symbol <- COAD_UCSC_Toil_tpm_dataset$COAD.RSEM.gene.tpm.symbol
gencode.v23.annotation <- COAD_UCSC_Toil_tpm_dataset$gencode.v23.annotation
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset$COAD.pheno
#2.read TCGA_COAD RSEM normal count Log2(x+1) data ####
COAD_RSEM_norm_count_Log2_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_Pan_TOIL_RSEM_norm_count/COAD_UCSC_RSEM_norm_count_Log2(x+1)_hugo_dataset.rds")
COAD.RSEM.gene.norm.count.Log2 <- COAD_RSEM_norm_count_Log2_dataset$`COAD.RSEM.gene.norm_count_Log2(x+1)`
COAD.pheno <- COAD_RSEM_norm_count_Log2_dataset$COAD.pheno










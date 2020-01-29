#### 2. Extract TCGA COAD data
### 1.Find subsset tcga data 
table(clinical_molecular_public_all$dataset)
clinical_molecular_public_tcga <- clinical_molecular_public_all[clinical_molecular_public_all$dataset == "tcga",]
table(clinical_molecular_public_tcga$msi)
### add rownames as the sample.tcga
sample.tcga <- as.character(clinical_molecular_public_tcga$sample)
rownames(clinical_molecular_public_tcga) <- gsub("-",".",sample.tcga)
head(clinical_molecular_public_tcga)
# add sampleID for table merge
clinical_molecular_public_tcga$sampleID <- gsub("-",".",sample.tcga)
### 2.Import TCGA COAD dataset clinicaldata 
COAD_UCSC_Toil_tpm_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset.rds")
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset$COAD.pheno
sample.COAD <- rownames(COAD.pheno)
sample.COAD.clean <- gsub('.{3}$', '', sample.COAD)
# add sampleID for table merge
COAD.pheno$sampleID <- sample.COAD.clean
# how many samples had the annotation 
table(COAD.pheno$sampleID %in% clinical_molecular_public_tcga$sampleID)
### 3.Merge two pheno table
COAD.pheno.merge.syn2623706 <- dplyr::left_join(COAD.pheno, clinical_molecular_public_tcga, by = "sampleID")
rownames(COAD.pheno.merge.syn2623706) <- rownames(COAD.pheno)
# add this mew pheno to the dataset
COAD_UCSC_Toil_tpm_dataset$COAD.pheno.merge.syn2623706 <- COAD.pheno.merge.syn2623706
saveRDS(COAD_UCSC_Toil_tpm_dataset,file = "/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
# test
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706 <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno.merge.syn2623706
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno

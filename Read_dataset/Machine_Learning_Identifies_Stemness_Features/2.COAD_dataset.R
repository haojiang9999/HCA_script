#### Extract COAD data
### 1.TCGA.Kallisto.cibersort_data.R
table(TCGA.Kallisto.cibersort.relative$CancerType)
COAD.cibersort.relative <- TCGA.Kallisto.cibersort.relative[TCGA.Kallisto.cibersort.relative$CancerType == "COAD",]
sum(COAD.cibersort.relative$B.cells.naive)
#3:24 cols
### For cols 3 -> 24 All results are reported as relative fractions normalized to 1 across all cell subsets
rowSums(COAD.cibersort.relative[,3:24])
### For "TotalLeukocyte" was  ESTIMATE applied to DNA methylation data
#We also obtained absolute estimates by scaling their relative abundance by overall
#leukocyte infiltration in each tumor as determined by ESTIMATE applied to DNA methylation data 
summary(COAD.cibersort.relative$TotalLeukocyte)
### 2.Purity_Ploidy using ESTIMATE
# Paper :Inferring tumour purity and stromal and immune cell admixture from expression data
# Modify the sample names
Purity_Ploidy_sampleID <- as.character(Purity_Ploidy_All_Samples_9_28_16$sample)
sampleID <- substr(Purity_Ploidy_sampleID,1,nchar(Purity_Ploidy_sampleID)-12)
Purity_Ploidy_All_Samples_9_28_16$SampleID <- gsub("-",".",sampleID)
### Check COAD data number
table(COAD.cibersort.relative$SampleID %in% Purity_Ploidy_All_Samples_9_28_16$SampleID)
table(Purity_Ploidy_All_Samples_9_28_16$SampleID %in% COAD.cibersort.relative$SampleID)
### 3.Merge two tables
COAD.cibersort.Purity <- dplyr::left_join(COAD.cibersort.relative,Purity_Ploidy_All_Samples_9_28_16, by = "SampleID")

COAD.cibersort.Purity$SampleID
cibersort.PurityID <- substr(COAD.cibersort.Purity$SampleID,1,nchar(COAD.cibersort.Purity$SampleID)-1)
### 4.Load COAD dataset
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706 <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno.merge.syn2623706
COAD.pheno$rownames <- rownames(COAD.pheno)
table(rownames(COAD.pheno) %in% cibersort.PurityID)
rownames(COAD.pheno) %in% dupID
dupID %in% rownames(COAD.pheno) 
dupID <- cibersort.PurityID[duplicated(cibersort.PurityID)]
COAD.cibersort.Purity[duplicated(cibersort.PurityID),]
COAD.cibersort.Purity[grep("TCGA.A6.6781.01",COAD.cibersort.Purity$SampleID),]
### 5.Merge the data to COAD.pheno
#### check every sample?
#### I have to remove duplicated samples by hand 
# I keep ones had more information or just Vial = "A"
#dupID
View(COAD.cibersort.Purity[grep(paste(dupID,collapse="|"),COAD.cibersort.Purity$SampleID),])
COAD.cibersort.Purity.unique <- COAD.cibersort.Purity[-c(12,14,19,32,41,43,47,50,54,59,72,78,80,373),]
# rownames same to COAD.pheno
rownames(COAD.cibersort.Purity.unique)<- substr(COAD.cibersort.Purity.unique$SampleID,1,nchar(COAD.cibersort.Purity.unique$SampleID)-1)
COAD.cibersort.Purity.unique$rownames <- rownames(COAD.cibersort.Purity.unique)
COAD.pheno.merge.all1 <- dplyr::left_join(COAD.pheno,COAD.cibersort.Purity.unique, by = "rownames")
### 6.Save the data to COAD dataset
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno.merge.all1 <- COAD.pheno.merge.all1
saveRDS(COAD_UCSC_Toil_tpm_dataset_merge.syn2623706, file = "/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")

COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno.merge.all1

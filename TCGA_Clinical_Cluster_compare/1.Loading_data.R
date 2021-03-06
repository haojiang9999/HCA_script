#### 1.Load data
### 1)Loading clinical data
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706 <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
COAD.pheno.merge.all1 <- COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno.merge.all1
### 2)Loading cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
### 3)Merge COAD.pheno.merge.all1 and Cluster.df
Plot.df <- dplyr::left_join(Cluster.df, COAD.pheno.merge.all1, by = "rownames")
### 4)Merge new clinical and molecular data
# read more clinical and molecular data
TCGA.COAD.ClinicalV2 <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_COAD_clinical_data/2020_2_10_TCGA.COAD.ClinicalV2.rds")
Plot.df <- dplyr::left_join(Plot.df, TCGA.COAD.ClinicalV2, by = "rownames")

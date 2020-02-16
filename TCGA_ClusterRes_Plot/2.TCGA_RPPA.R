#### 2.TCGA_RPPA.R
# reverse-phase protein array (RPPA)
## 1.Read dataset
COAD_TCGA_RPPA_clean_xena_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/RPPA_n_7744/COAD_TCGA_RPPA_clean_xena_dataset.rds")
COAD.RPPA.pancan.clean.xena <- COAD_TCGA_RPPA_clean_xena_dataset$COAD.RPPA.pancan.clean.xena
COAD.RPPA.pancan.clean.xena[1:5,1:5]

#### Read scRNA colon cancer
### 1.Tang colon cancer
Tang_scColon_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Tang_scColon_cancer_dataset/Tang_scColon_dataset.rds")
Tang.colon.cancer.FPKM.500 <- Tang_scColon_dataset$FPKM.500.exp
Tang.colon.cancer.TPM.688 <- Tang_scColon_dataset$TPM.688.exp
Tang.colon.cancer.FPKM.500.pheno <- Tang_scColon_dataset$FPKM.500.pheno
Tang.colon.cancer.TPM.688.pheno <- Tang_scColon_dataset$TPM.688.pheno
### 2.RCA colon cancer
RCA_epithelial_COUNT_normal_uniGene_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_epithelial_COUNT_normal_uniGene_dataset.rds")
NM.epithelial.COUNT.normal.uniGene <- RCA_epithelial_COUNT_normal_uniGene_dataset$NM.epithelial.COUNT.normal.uniGene
Tumor.epithelial.COUNT.normal.uniGene <- RCA_epithelial_COUNT_normal_uniGene_dataset$Tumor.epithelial.COUNT.normal.uniGene
NM.Epi.phenoType <- RCA_epithelial_COUNT_normal_uniGene_dataset$NM.Epi.phenoType
Tumor.Epi.phenoType <- RCA_epithelial_COUNT_normal_uniGene_dataset$tumor.Epi.phenoType

## convert dgCMatrix to dataframe
#NM.epithelial.COUNT.normal.uniGene <- as.data.frame(as.matrix(NM.epithelial.COUNT.normal.uniGene)) 
#Tumor.epithelial.COUNT.normal.uniGene <- as.data.frame(as.matrix(Tumor.epithelial.COUNT.normal.uniGene)) 

### 3.Colon organoid
Colon_Orgnoid_uniGene_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/3D_organoid/3D_Organoid/Colon_Orgnoid_uniGene_dataset.rds")
Orgnoid_expression.uniGene.TPM <- Colon_Orgnoid_uniGene_dataset$Orgnoid_expression.uniGene.TPM
Orgnoid.colon.AUC.IC50 <- Colon_Orgnoid_uniGene_dataset$Drug_AUC_IC50




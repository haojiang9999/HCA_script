## read the RCA FPKM epithelial dataset
filePath <- "/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_epithelial_FPKM_dataset.rds"
RCA_epithelial_FPKM_dataset.rds <- readRDS(filePath)
NM.epithelial.fpkm.exp <- RCA_epithelial_FPKM_dataset.rds$NM.epithelial.fpkm.exp
tumor.epithelial.fpkm.exp <- RCA_epithelial_FPKM_dataset.rds$tumor.epithelial.fpkm.exp
NM.Epi.phenoType <- RCA_epithelial_FPKM_dataset.rds$NM.Epi.phenoType
tumor.Epi.phenoType <- RCA_epithelial_FPKM_dataset.rds$tumor.Epi.phenoType
GeneAnno <- RCA_epithelial_FPKM_dataset.rds$GeneAnno

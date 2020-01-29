## read the RCA COUNT epithelial dataset
filePath <- "/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_epithelial_COUNT_dataset.rds"
RCA_epithelial_COUNT_dataset.rds <- readRDS(filePath)
NM.epithelial.COUNT.exp <- RCA_epithelial_COUNT_dataset.rds$NM.epithelial.COUNT.exp
tumor.epithelial.COUNT.exp <- RCA_epithelial_COUNT_dataset.rds$tumor.epithelial.COUNT.exp
NM.Epi.phenoType <- RCA_epithelial_COUNT_dataset.rds$NM.Epi.phenoType
tumor.Epi.phenoType <- RCA_epithelial_COUNT_dataset.rds$tumor.Epi.phenoType
GeneAnno <- RCA_epithelial_COUNT_dataset.rds$GeneAnno

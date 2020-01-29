#### Read gene name filtered lncRNA
filePath <- "/stor/jianghao/GEO/Cancer/GSE81861_human_colorectal_tumors/R_GSE81861_human_colorectal_tumors/epi_stem_DE_TvN.lnc.NameFiltered_fc_0.2.csv"
epi_stem_DE_TvN.lnc.NameFiltered_fc_0.2 <- read.csv(filePath)  
# Read gene expression dataset 
RCA_epithelial_FPKM_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_epithelial_FPKM_dataset.rds")
NM.epithelial.fpkm.exp <- RCA_epithelial_FPKM_dataset$NM.epithelial.fpkm.exp
tumor.epithelial.fpkm.exp <- RCA_epithelial_FPKM_dataset$tumor.epithelial.fpkm.exp
NM.Epi.phenoType <- RCA_epithelial_FPKM_dataset$NM.Epi.phenoType
tumor.Epi.phenoType <- RCA_epithelial_FPKM_dataset$tumor.Epi.phenoType
GeneAnno <- RCA_epithelial_FPKM_dataset$GeneAnno
## stem cell compare
samNM.stemTA <- rownames(NM.Epi.phenoType)[NM.Epi.phenoType$Epi_cellTypes == "stemTA"]
samTumor.stemTA <- rownames(tumor.Epi.phenoType)[tumor.Epi.phenoType$Epi_cellTypes == "stemTA"]
## Read pQ normalized data
RCA_epithelial.fpkm.exp.stemTA.pQ <- readRDS("/data8t_4/JH/MyJobs/ColonCancer_lncRNA/RCA_epithelial_analysis/RCA_epithelial.fpkm.exp.stemTA.pQ.rds")
RCA_all.fpkm.exp.pQ <- readRDS("/data8t_4/JH/MyJobs/ColonCancer_lncRNA/RCA_epithelial_analysis/RCA_all.fpkm.exp.pQ.rds")

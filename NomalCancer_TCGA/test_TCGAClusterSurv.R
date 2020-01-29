## test TCGAClusterSurv
Cor.tumor <- Cor.Res.CV2000$Cor.merged
Cor.tumor<- Cor.tumorCor.tumor.cv10
Input.tb <- Cor.tumor
hclust.res <- hcTumor.cv10
COAD.pheno
Col.anno <- COAD.pheno[TumorID,]
colnames(COAD.pheno)
source("/data8t_4/JH/MyJobs/1_R_script/TCGA_plot/TCGAClusterSurv.R")
TCGAClusterSurv(Input.tb = Input.tb, hclust.res = hclust.res, Col.anno = Col.anno, k = 3,
                col_anno = c("sampleTypes","histological_type"))

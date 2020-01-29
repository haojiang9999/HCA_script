#### Top varience features
TumorID <- rownames(COAD.pheno[COAD.pheno$sampleTypes == "Tumor",])
Pheno.merged.tumor <- COAD.pheno[TumorID,c("sampleTypes","histological_type")]
## cv2500 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/JH_SC3_cluster.R")
Cor.tumor <- Cor.Res.CV2500$Cor.merged
Cor.tumor<- Cor.tumor[,TumorID]
###
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R") 
Cor.tumor.top30<- TopCV(Cor.tumor, TopN = 30, MARGIN = 1)
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
test1<- heatmap.JH(Cor.tumor.top30,show_colnames = F,cluster_rows = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = T)
TCGAClusterSurv(Input.tb = Cor.tumor.top30, hclust.res = test1$tree_col, Col.anno = Col.anno, k = 2)

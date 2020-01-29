#### Function for TCGA Cluster survival analysis
TCGAClusterSurv <- function(Input.tb, hclust.res, Col.anno, k = 2){
  # Input.tb
  #Cor.tumor <- Cor.Res.CV2000$Cor.merged
  #Cor.tumor<- Cor.tumor[,TumorID]
  #Input.tb <- Cor.tumor
  #hclust.res <- hcTumor.2000
  #Col.anno <- COAD.pheno[TumorID,c("sampleTypes","histological_type","OS","OS.time")]
  ### 1.cut the tree
  cutree.res <- cutree(hclust.res, k = k)
  Col.anno.cluster <- cbind(Col.anno[,c("sampleTypes","histological_type")],as.character(cutree.res))
  ### 2.Cut tree cluster result
  source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
  plot(as.dendrogram(hclust.res)) 
  heatmap.JH(Input.tb,show_colnames = F,cluster_cols = hclust.res,
             annotation_col =Col.anno.cluster)
  ### 3.Build the survival table
  sampleID <- names(cutree.res)
  surTime <- Col.anno[sampleID,c("OS","OS.time")]
  sur.df <- cbind(cutree.res,surTime)
  ### 4.Survival analysis
  if (!require("survminer")) 
    BiocManager::install("survminer")
  require(survminer)
  require("survival")
  fit <- survfit(Surv(OS.time, OS) ~ cutree.res, data = sur.df)
  # Drawing curves
  ggsurvplot(fit,pval = T,pval.method = T)
}

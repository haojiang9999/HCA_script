TCGAClusterSurv.2<- function(Input.tb, hclust.res, group.res, Col.anno, col_anno = c("sampleTypes","histological_type")){
  # Input.tb
  #Cor.tumor <- Cor.Res.CV2000$Cor.merged
  #Cor.tumor<- Cor.tumor[,TumorID]
  #Input.tb <- Cor.tumor
  #hclust.res <- hcTumor.2000
  #Col.anno <- COAD.pheno[TumorID,c("sampleTypes","histological_type","OS","OS.time")]
  ### 1.cut the tree
  #cutree.res <- cutree(hclust.res, k = k)
  Col.anno.cluster <- cbind(Col.anno[,col_anno],as.character(group.res))
  ### 2.Cut tree cluster result
  source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
  plot(as.dendrogram(hclust.res)) 
  heatmap.JH(Input.tb,show_colnames = F,cluster_cols = hclust.res,cluster_rows = F,
             annotation_col =Col.anno.cluster)
  ### 3.Build the survival table
  sampleID <- colnames(Input.tb)
  surTime <- Col.anno[sampleID,c("OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")]
  sur.df <<- cbind(group.res,surTime)
  ### 4.Survival analysis
  if (!require("survminer")) 
    BiocManager::install("survminer")
  require(survminer)
  require("survival")
  fitOS <- survfit(Surv(OS.time, OS) ~ group.res, data = sur.df)
  fitDSS <- survfit(Surv(DSS.time, DSS) ~ group.res, data = sur.df)
  fitDFI <- survfit(Surv(DFI.time, DFI) ~ group.res, data = sur.df)
  fitPFI <- survfit(Surv(PFI.time, PFI) ~ group.res, data = sur.df)
  # Drawing curves
  #ggsurvplot(fitOS,pval = T,pval.method = T,title = "Over all survival(OS)")
  #ggsurvplot(fitDSS,pval = T,pval.method = T,title = "Disease-specific survival (DSS)")
  #ggsurvplot(fitDFI,pval = T,pval.method = T,main = "disease-free interval (DFI)")
  #ggsurvplot(fitDFI,pval = T,pval.method = T,main = "Progression-free interval (PFI) ")
  title <- c("Over all survival(OS)", "Disease-specific survival (DSS)",
             "disease-free interval (DFI)","Progression-free interval (PFI)")
  ggsurvplot_list(list(fitOS,fitDSS,fitDFI,fitPFI), sur.df, title = title, risk.table=T)
}

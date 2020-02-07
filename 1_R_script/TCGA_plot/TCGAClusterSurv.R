#### Function for TCGA Cluster survival analysis
TCGAClusterSurv <- function(Input.tb, hclust.res, Col.anno, k = 2, col_anno = c("sampleTypes","histological_type")){
  # Input.tb
  #Cor.tumor <- Cor.Res.CV2000$Cor.merged
  #Cor.tumor<- Cor.tumor[,TumorID]
  #Input.tb <- Cor.tumor
  #hclust.res <- hcTumor.2500
  #Col.anno <- COAD.pheno[TumorID,c("sampleTypes","histological_type","OS","OS.time")]
  #Col.anno <- COAD.pheno[TumorID,]
  ### 1.cut the tree
  cutree.res <- cutree(hclust.res, k = k)
  if (!require(WGCNA)){
    source("http://bioconductor.org/biocLite.R");
    biocLite(c("impute", "GO.db", "preprocessCore")); 
    install.packages("WGCNA");
  }
  require(WGCNA)  
  dynamicColors = labels2colors(cutree.res) # convert label to color
  Col.anno.cluster <- cbind(Col.anno[,col_anno],cluster.Res = dynamicColors)
  cluster.Res.color <- unique(dynamicColors)
  names(cluster.Res.color) <- cluster.Res.color
  mycolors <- list(cluster.Res = cluster.Res.color)
  ### 2.Cut tree cluster result
  source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
  plot(as.dendrogram(hclust.res)) 
  heatmap.JH(Input.tb,show_colnames = F,cluster_cols = hclust.res,cluster_rows = F,
             annotation_col =Col.anno.cluster, 
             annotation_colors = mycolors)
  ### 3.Build the survival table
  sampleID <- names(cutree.res)
  surTime <- Col.anno[sampleID,c("OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")]
  sur.df <<- cbind(dynamicColors,surTime)
  ### 4.Survival analysis
  if (!require("survminer")) 
    BiocManager::install("survminer")
  require(survminer)
  require("survival")
  fitOS <- survfit(Surv(OS.time, OS) ~ dynamicColors, data = sur.df)
  fitDSS <- survfit(Surv(DSS.time, DSS) ~ dynamicColors, data = sur.df)
  fitDFI <- survfit(Surv(DFI.time, DFI) ~ dynamicColors, data = sur.df)
  fitPFI <- survfit(Surv(PFI.time, PFI) ~ dynamicColors, data = sur.df)
  # Drawing curves
  ggsurvplot(fitOS,pval = T,pval.method = T,title = "Over all survival(OS)",  palette = sort(as.character(cluster.Res.color)))
  #ggsurvplot(fitDSS,pval = T,pval.method = T,title = "Disease-specific survival (DSS)")
  #ggsurvplot(fitDFI,pval = T,pval.method = T,main = "disease-free interval (DFI)")
  #ggsurvplot(fitDFI,pval = T,pval.method = T,main = "Progression-free interval (PFI) ")
  title <- c("Over all survival(OS)", "Disease-specific survival (DSS)",
             "disease-free interval (DFI)","Progression-free interval (PFI)")
  ggsurvplot_list(list(fitOS,fitDSS,fitDFI,fitPFI), sur.df, title = title, risk.table=T,
                  palette = sort(as.character(cluster.Res.color)))
}

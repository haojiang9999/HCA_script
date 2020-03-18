#### 3.For_loop_silhouette_score_calculation.R
TopNs <- seq(1000,6000,by = 100)
si.sum <- list()
for(i in 1:length(TopNs)){
  ### 2.Find genes had high coefficience-varience
  # i=4
  TopN = TopNs[i] ###########################################
  source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
  scReference.list.TopN <- lapply(scReference.list.V8.GI, function(x){
    TopCV(x, TopN = TopN, MARGIN = 1)
  })
  
  ### 3.log10(x+1) transformed scReference.list.TopN
  log10.scReference.list.TopN <- lapply(scReference.list.TopN, function(x){
    log10.x <- log10(x+1)
    return(log10.x)
  })
  
  
  ### 4.Correlation calculation 
  source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
  Cor.Res.TopN <- refCorMerge(COAD_tpm.Log10.expList, log10.scReference.list.TopN)
  Cor.merged.TopN <- Cor.Res.TopN$Cor.merged
  COAD.pheno <- COAD.pheno
  #Cor.merged.TopN
  
  ### 5.Correlation normalization steps
  ## 1)The first Normalization
  #normal the correlation across samples
  Trans.Rang1.TopN<- base::apply(Cor.merged.TopN, 1, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  Trans.Rang1.TopN <- t(Trans.Rang1.TopN)
  
  ## 2)The second Normalization
  # Normalize the correlation across each cell type cluster; It's the weight of each cell type cluster
  cellTypesName <- rownames(Trans.Rang1.TopN)
  #grep("Adult", cellTypesName)
  Trans.Rang1.TopN.Normalized <- Trans.Rang1.TopN
  Trans.Rang1.TopN.Normalized[grep("Adult", cellTypesName),]<- base::apply(Trans.Rang1.TopN.Normalized[grep("Adult", cellTypesName),], 2, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  Trans.Rang1.TopN.Normalized[grep("Fetal", cellTypesName),]<- base::apply(Trans.Rang1.TopN.Normalized[grep("Fetal", cellTypesName),], 2, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  
  ### 6.Distance calculation
  Trans.Rang1.TopN.Normalized
  TopN.Normalized.dist <- dist(t(Trans.Rang1.TopN.Normalized), method="euclidean")
  
  ### 7.Calculate the silhouette
  library("cluster")
  #cbind(rownames(COAD.pheno),colnames(Trans.Rang1.TopN.Normalized))
  NTcluster <- COAD.pheno$sampleTypes
  NTcluster <- as.numeric(NTcluster)
  #table(NTcluster)
  names(NTcluster) <- rownames(COAD.pheno)
  # Calculate silhouette value for each sample
  sil <- silhouette(NTcluster, TopN.Normalized.dist)
  # Summary of silhouette analysis
  # Average silhouette width of each cluster
  si.sum[[i]] <- summary(sil)$clus.avg.widths
  library(factoextra)
  print(paste("Top",TopNs[i]))
  p <- fviz_silhouette(sil,title=paste("Top",TopNs[i]))
  
  print(p)

}

## Merge table
si.sum.merge <- do.call(rbind,si.sum)
si.sum.merge <- as.data.frame(si.sum.merge)
si.sum.merge$TopNCV <- TopNs 
## ggplot
library(ggplot2)
library(reshape2)
si.sum.merge.m <- melt(si.sum.merge,id.vars = "TopNCV")
ggplot(si.sum.merge.m, aes(TopNCV,value,group = variable)) + 
  geom_line() +
  geom_point(aes(color=variable)) +
  geom_vline(xintercept = 2500,linetype="dashed")

### 500 interval
si.sum.merge.500 <- si.sum.merge[si.sum.merge$TopNCV %in% seq(1000,6000,by = 500),]
si.sum.merge.500.m <- melt(si.sum.merge.500,id.vars = "TopNCV")
ggplot(si.sum.merge.500.m, aes(TopNCV,value,group = variable)) + 
  geom_line() +
  geom_point(aes(color=variable)) +
  geom_vline(xintercept = 2500,linetype="dashed")

## Clustering Validation Statistics
#http://www.sthda.com/english/wiki/wiki.php?id_contents=7952#silhouette-plot-for-hierarchical-clustering
sil
sil$silinfo
library(factoextra)
fviz_silhouette(sil)
# Summary of silhouette analysis
si.sum <- summary(sil)
# Average silhouette width of each cluster
si.sum$clus.avg.widths
# The total average (mean of all individual silhouette widths)
si.sum$avg.width
# The size of each clusters
si.sum$clus.sizes
x <- fviz_silhouette(sil)
# Summary of silhouette analysis
si.sum <- summary(sil)
# Average silhouette width of each cluster
si.sum$clus.avg.widths

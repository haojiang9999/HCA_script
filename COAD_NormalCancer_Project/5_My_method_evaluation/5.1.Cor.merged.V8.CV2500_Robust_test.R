#### 5.1Cor.merged.V8.CV2500_Robust_test.R
### 1.Loading data
Cor.merged.V8.CV2500.dataset <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/5_My_method_evaluation/Cor.merged.V8.CV2500.dataset.rds")
Cor.merged.V8.CV2500 <- Cor.merged.V8.CV2500.dataset$Cor.merged.V8.CV2500
COAD.pheno <- Cor.merged.V8.CV2500.dataset$COAD.pheno
NTcluster <- COAD.pheno$sampleTypes
NTcluster <- as.numeric(NTcluster)
#table(NTcluster)
names(NTcluster) <- rownames(COAD.pheno)
### 2.Remove scReference corelation
dim(Cor.merged.V8.CV2500)
## remove 1 item
si.average.merge.1 <- list()
for(i in 1:39){
  # i=1
  Cor.merged.V8.CV2500_sub <- Cor.merged.V8.CV2500[-i,]
  Cor.merged.V8.CV2500_sub.dist <- dist(t(Cor.merged.V8.CV2500_sub), method="euclidean")
  require("cluster")
  sil.merge <- silhouette(NTcluster, Cor.merged.V8.CV2500_sub.dist)
  si.average.merge.1[[i]] <- summary(sil.merge)$clus.avg.widths
}
## Merge table
si.average.merge.1.tb <- do.call(rbind,si.average.merge.1)
si.average.merge.1.tb <- as.data.frame(si.average.merge.1.tb)
si.average.merge.1.tb$RMnum <- 1:39 
## ggplot
library(ggplot2)
library(reshape2)
si.average.merge.1.tb.m <- melt(si.average.merge.1.tb,id.vars = "RMnum")
ggplot(si.average.merge.1.tb.m, aes(RMnum,value,group = variable)) + 
  geom_line() +
  geom_point(aes(color=variable))
library(pheatmap)
COAD.pheno_sub <- as.data.frame(COAD.pheno[,c("sampleTypes","msi")])
pheatmap::pheatmap(Cor.merged.V8.CV2500_sub, annotation_col = COAD.pheno_sub)

Cor.merged.V8.CV2500

library(factoextra)
fviz_silhouette(sil)


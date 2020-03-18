#### 5.2.Cor.merged.V8.CV2500.cor.Scaled2s_Robust_test.R
### 1.Loading data
Cor.merged.V8.CV2500.scaled.Twice.dataset <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/5_My_method_evaluation/Cor.merged.V8.CV2500.scaled.Twice.dataset.rds")
Trans.Rang1.CV2500.Normalized <- Cor.merged.V8.CV2500.scaled.Twice.dataset$Trans.Rang1.CV2500.Normalized
COAD.pheno <- Cor.merged.V8.CV2500.scaled.Twice.dataset$COAD.pheno
NTcluster <- COAD.pheno$sampleTypes
NTcluster <- as.numeric(NTcluster)

### 2.Remove scReference corelation
dim(Trans.Rang1.CV2500.Normalized)
## remove 1 item
si.average <- list()
for(i in 1:39){
  # i=1
  Trans.Rang1.CV2500.Normalized_sub <- Trans.Rang1.CV2500.Normalized[-i,]
  Trans.Rang1.CV2500.Normalized_sub.dist <- dist(t(Trans.Rang1.CV2500.Normalized_sub), method="euclidean")
  require("cluster")
  sil.merge <- silhouette(NTcluster, Trans.Rang1.CV2500.Normalized_sub.dist)
  si.average[[i]] <- summary(sil.merge)$clus.avg.widths
  COAD.pheno_sub <- as.data.frame(COAD.pheno[,c("sampleTypes","msi")])
  pheatmap::pheatmap(Trans.Rang1.CV2500.Normalized_sub, 
                     annotation_col = COAD.pheno_sub,
                     cluster_rows = F)
}
## Merge table
si.average.tb <- do.call(rbind,si.average)
si.average.tb <- as.data.frame(si.average.tb)
si.average.tb$RMnum <- 1:39 
## ggplot
library(ggplot2)
library(reshape2)
si.average.tb.m <- melt(si.average.tb,id.vars = "RMnum")
ggplot(si.average.tb.m, aes(RMnum,value,group = variable)) + 
  geom_line() +
  geom_point(aes(color=variable))


pheatmap::pheatmap(Trans.Rang1.CV2500.Normalized, 
                   annotation_col = COAD.pheno_sub,
                   cluster_rows = F)
Trans.Rang1.CV2500.Normalized.dist <- dist(t(Trans.Rang1.CV2500.Normalized), method="euclidean")
sil.merge <- silhouette(NTcluster, Trans.Rang1.CV2500.Normalized.dist)
summary(sil.merge)$clus.avg.widths
library(factoextra)
fviz_silhouette(sil.merge)

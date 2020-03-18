#### 5.3.Cor.merged.V8.CV2500.cor.Scaled2s_Robust_test.rm2.R
### 1.Loading data
Cor.merged.V8.CV2500.scaled.Twice.dataset <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/5_My_method_evaluation/Cor.merged.V8.CV2500.scaled.Twice.dataset.rds")
Trans.Rang1.CV2500.Normalized <- Cor.merged.V8.CV2500.scaled.Twice.dataset$Trans.Rang1.CV2500.Normalized
COAD.pheno <- Cor.merged.V8.CV2500.scaled.Twice.dataset$COAD.pheno
NTcluster <- COAD.pheno$sampleTypes
NTcluster <- as.numeric(NTcluster)
### 2.Remove scReference corelation
dim(Trans.Rang1.CV2500.Normalized)
cellTypesName <- rownames(Trans.Rang1.CV2500.Normalized)
## Combination of select 2 items from cellTypes
# prod(1:39)/(prod(1:2)*prod(1:37))
c2 <- t(combn(cellTypesName, 2))
class(c2)
## remove 2 items
si.average <- list()
for(i in 1:nrow(c2) ){
  # i=1
  rmCellTypes <- c2[i,]
  Trans.Rang1.CV2500.Normalized_sub <- Trans.Rang1.CV2500.Normalized[rownames(Trans.Rang1.CV2500.Normalized) != rmCellTypes,]
  Trans.Rang1.CV2500.Normalized_sub.dist <- dist(t(Trans.Rang1.CV2500.Normalized_sub), method="euclidean")
  require("cluster")
  sil.merge <- silhouette(NTcluster, Trans.Rang1.CV2500.Normalized_sub.dist)
  si.average[[i]] <- summary(sil.merge)$clus.avg.widths
  #COAD.pheno_sub <- as.data.frame(COAD.pheno[,c("sampleTypes","msi")])
  #pheatmap::pheatmap(Trans.Rang1.CV2500.Normalized_sub, 
  #                   annotation_col = COAD.pheno_sub,
  #                   cluster_rows = F)
}
si.average.tb <- do.call(rbind,si.average)
si.average.tb <- as.data.frame(si.average.tb)



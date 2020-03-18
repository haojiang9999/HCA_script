#### 4.1.CV2500_Random_remove_scReference_Prepare.R
### 1.Loading data
## 1)loading COAD_tpm_log10.expList ###
COAD_tpm.Log10.expList.dataset <- readRDS("COAD_tpm.Log10.expList.dataset.rds")
COAD_tpm.Log10.expList <- COAD_tpm.Log10.expList.dataset$COAD_tpm.Log10.expList
COAD.pheno <- COAD_tpm.Log10.expList.dataset$COAD.pheno

### 2.Top CV2500 scReference panel building
## 1)Read the Step1 build scReference table
# Only my cluster using reference
Tang.Adult.colon.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Adult.colon.ref.rds")
Tang.Fetal.GI.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Fetal.GI.ref.rds")
## 2)Build_scReference_list
# Remove Adult_Immune immuno cell from Tang.Adult.colon.ref
#colnames(Tang.Adult.colon.ref)
scReference.list.V8.GI <- list(Tang.Adult.colon = Tang.Adult.colon.ref[,-which(names(Tang.Adult.colon.ref) %in% "Adult_Immune")],
                               Tang.Fetal.GI = Tang.Fetal.GI.ref)


source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scReference.list.CV2500 <- lapply(scReference.list.V8.GI, function(x){
  TopCV(x, TopN = 2500, MARGIN = 1)
})

### 3.log10(x+1) transformed scReference.list.TopN
log10.scReference.list.CV2500 <- lapply(scReference.list.CV2500, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})


### 4.Correlation calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
Cor.Res.CV2500 <- refCorMerge(COAD_tpm.Log10.expList, log10.scReference.list.CV2500)
Cor.merged.CV2500 <- Cor.Res.CV2500$Cor.merged
COAD.pheno <- COAD.pheno
#Cor.merged.TopN

### 5.Correlation normalization steps
## 1)The first Normalization
#normal the correlation across samples
Trans.Rang1.CV2500<- base::apply(Cor.merged.CV2500, 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.CV2500 <- t(Trans.Rang1.CV2500)

## 2)The second Normalization
# Normalize the correlation across each cell type cluster; It's the weight of each cell type cluster
cellTypesName <- rownames(Trans.Rang1.CV2500)
#grep("Adult", cellTypesName)
Trans.Rang1.CV2500.Normalized <- Trans.Rang1.CV2500
Trans.Rang1.CV2500.Normalized[grep("Adult", cellTypesName),]<- base::apply(Trans.Rang1.CV2500.Normalized[grep("Adult", cellTypesName),], 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.CV2500.Normalized[grep("Fetal", cellTypesName),]<- base::apply(Trans.Rang1.CV2500.Normalized[grep("Fetal", cellTypesName),], 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
## 3) save result
saveRDS(Trans.Rang1.CV2500.Normalized, file = "Cor.merged.V8.CV2500.scaled.Twice.dataset.rds")

### 6.Distance calculation
#Trans.Rang1.CV2500.Normalized
Trans.Rang1.CV2500.Normalized.dist <- dist(t(Trans.Rang1.CV2500.Normalized), method="euclidean")

### 7.Calculate the silhouette
library("cluster")
#cbind(rownames(COAD.pheno),colnames(Trans.Rang1.TopN.Normalized))
NTcluster <- COAD.pheno$sampleTypes
NTcluster <- as.numeric(NTcluster)
#table(NTcluster)
names(NTcluster) <- rownames(COAD.pheno)
# Calculate silhouette value for each sample
sil <- silhouette(NTcluster, Trans.Rang1.CV2500.Normalized.dist)
# Summary of silhouette analysis
# Average silhouette width of each cluster
si.sum <- summary(sil)$clus.avg.widths
library(factoextra)
#print(paste("Top",TopNs[i]))
p <- fviz_silhouette(sil,title=paste("Top",TopNs[i]))

print(p)


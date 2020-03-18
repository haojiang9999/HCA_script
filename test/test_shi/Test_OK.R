## Test_OK.R

###loading
## 1)loading COAD_tpm_log10.expList ###
COAD_tpm.Log10.expList.dataset <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/5_My_method_evaluation/COAD_tpm.Log10.expList.dataset.rds")
COAD_tpm.Log10.expList <- COAD_tpm.Log10.expList.dataset$COAD_tpm.Log10.expList
COAD.pheno <- COAD_tpm.Log10.expList.dataset$COAD.pheno
Tang.Adult.colon.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Adult.colon.ref.rds")
Tang.Fetal.GI.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Fetal.GI.ref.rds")
## 2)Build_scReference_list
# Remove Adult_Immune immuno cell from Tang.Adult.colon.ref
#colnames(Tang.Adult.colon.ref)
scReference.list.V8.GI <- list(Tang.Adult.colon = Tang.Adult.colon.ref[,-which(names(Tang.Adult.colon.ref) %in% "Adult_Immune")],
                               Tang.Fetal.GI = Tang.Fetal.GI.ref)


"/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/5_My_method_evaluation/"
### 2.Find genes had high coefficience-varience
# i=4

source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
scReference.list.TopN <- lapply(scReference.list.V8.GI, function(x){
  TopCV(x, TopN = 2500, MARGIN = 1)
})

### 3.log10(x+1) transformed scReference.list.TopN
log10.scReference.list.TopN <- lapply(scReference.list.TopN, function(x){
  log10.x <- log10(x+1)
  return(log10.x)
})


### 4.Correlation calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
Cor.Res.TopN <- refCorMerge(COAD_tpm.Log10.expList, scReference.list.TopN)
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
si.sum <- summary(sil)$clus.avg.widths
library(factoextra)
#print(paste("Top",TopNs[i]))
p <- fviz_silhouette(sil)

print(p)
COAD.pheno_sub <- as.data.frame(COAD.pheno[,c("sampleTypes","msi")])
pheatmap::pheatmap(Trans.Rang1.TopN.Normalized, annotation_col = COAD.pheno_sub)
pheatmap::pheatmap(Trans.Rang1.TopN.Normalized, )
#### 1.CV_value_choosing_test.R
### 1.Loading data
### scReference loading ###
## 1)Read the Step1 build scReference table
# Only my cluster using reference
Tang.Adult.colon.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Adult.colon.ref.rds")
Tang.Fetal.GI.ref <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step1_Normal_cell_reference_panel_construction/Tang.Fetal.GI.ref.rds")
## 2)Build_scReference_list
# Remove Adult_Immune immuno cell from Tang.Adult.colon.ref
#colnames(Tang.Adult.colon.ref)
scReference.list.V8.GI <- list(Tang.Adult.colon = Tang.Adult.colon.ref[,-which(names(Tang.Adult.colon.ref) %in% "Adult_Immune")],
                            Tang.Fetal.GI = Tang.Fetal.GI.ref)
## 3)loading COAD_tpm_log10.expList ###
COAD_tpm.Log10.expList.dataset <- readRDS("COAD_tpm.Log10.expList.dataset.rds")
COAD_tpm.Log10.expList <- COAD_tpm.Log10.expList.dataset$COAD_tpm.Log10.expList
COAD.pheno <- COAD_tpm.Log10.expList.dataset$COAD.pheno
#class(COAD_tpm.Log10.expList)

### 2.Find genes had high coefficience-varience
TopN = 1000 ###########################################
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

### 7.Calculate the Dunn
library(clValid)
#cbind(rownames(COAD.pheno),colnames(Trans.Rang1.TopN.Normalized))
NTcluster <- COAD.pheno$sampleTypes
NTcluster <- as.numeric(NTcluster)
#table(NTcluster)
names(NTcluster) <- rownames(COAD.pheno)
dunn.score <- dunn(TopN.Normalized.dist,NTcluster )

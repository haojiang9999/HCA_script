#### 1.1_Loading_2020_2_1_enviroment_V7_Tumor
# CutTree1.Loading_data.R
#### 1.Load dataset
###### read TCGA_COAD data ####
COAD_UCSC_Toil_tpm_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
COAD_RSEM_gene_tpm <- COAD_UCSC_Toil_tpm_dataset$COAD.RSEM.gene.tpm
gencode.v23.annotation <- COAD_UCSC_Toil_tpm_dataset$gencode.v23.annotation
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset$COAD.pheno.merge.syn2623706
#head(gencode.v23.annotation)
## convert ensembleID to symbol 
geneMatch <- match(rownames(COAD_RSEM_gene_tpm),gencode.v23.annotation$V1)
geneSymbol <- as.character(gencode.v23.annotation[geneMatch,]$V2)
COAD_tpm_symbol <- COAD_RSEM_gene_tpm
rownames(COAD_tpm_symbol) <- geneSymbol
#head(COAD_tpm_symbol)
## Read log10(x+1) transformed scReference 

#### Load scReference V5 
scReference.log10.CV <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step2_Merge_and_Filter_the_scReference/2020_2_1_scReferenceV7.log10.CV.rds")
scReference.list.log10.CV8000 <- scReference.log10.CV$scReference.list.log10.CV8000
scReference.list.log10.CV4500 <- scReference.log10.CV$scReference.list.log10.CV4500
scReference.list.log10.CV4000 <- scReference.log10.CV$scReference.list.log10.CV4000
scReference.list.log10.CV3500 <- scReference.log10.CV$scReference.list.log10.CV3500
scReference.list.log10.CV3000 <- scReference.log10.CV$scReference.list.log10.CV3000
scReference.list.log10.CV2500 <- scReference.log10.CV$scReference.list.log10.CV2500
scReference.list.log10.CV2000 <- scReference.log10.CV$scReference.list.log10.CV2000
#### 2.Distance calculation
### Distance calculation
##### Distance calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
##### Transform COAD data
#summary(colSums(COAD_tpm_symbol)) ### So its not log2 transformed
# Expression data transformation Log(x+1)
Log10.expList <- list(COAD_tpm_symbol = log10(COAD_tpm_symbol+1))
Cor.Res.CV8000 <- refCorMerge(Log10.expList, scReference.list.log10.CV8000)
Cor.Res.CV4500 <- refCorMerge(Log10.expList, scReference.list.log10.CV4500)
Cor.Res.CV4000 <- refCorMerge(Log10.expList, scReference.list.log10.CV4000)
Cor.Res.CV3500 <- refCorMerge(Log10.expList, scReference.list.log10.CV3500)
Cor.Res.CV3000 <- refCorMerge(Log10.expList, scReference.list.log10.CV3000)
Cor.Res.CV2500 <- refCorMerge(Log10.expList, scReference.list.log10.CV2500)
Cor.Res.CV2000 <- refCorMerge(Log10.expList, scReference.list.log10.CV2000)
### Find tumor samples
### only tumor samples
TumorID <- rownames(COAD.pheno[COAD.pheno$sampleTypes == "Tumor",])
COAD.pheno.tumor <- COAD.pheno[TumorID,c("sampleTypes","msi","cimp")]

#### 3.Normalization Tumor samples
### normalize across the features
Trans.Rang1.cv8000<- base::apply(Cor.Res.CV8000$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv4500<- base::apply(Cor.Res.CV4500$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv4000<- base::apply(Cor.Res.CV4000$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv3500<- base::apply(Cor.Res.CV3500$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv3000<- base::apply(Cor.Res.CV3000$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv2500<- base::apply(Cor.Res.CV2500$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv2000<- base::apply(Cor.Res.CV2000$Cor.merged[,TumorID], 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv8000 <- t(Trans.Rang1.cv8000)
Trans.Rang1.cv4500 <- t(Trans.Rang1.cv4500)
Trans.Rang1.cv4000 <- t(Trans.Rang1.cv4000)
Trans.Rang1.cv3500 <- t(Trans.Rang1.cv3500)
Trans.Rang1.cv3000 <- t(Trans.Rang1.cv3000)
Trans.Rang1.cv2500 <- t(Trans.Rang1.cv2500)
Trans.Rang1.cv2000 <- t(Trans.Rang1.cv2000)


# 1.Subset the Adult and fetal scReference
Trans.Rang1.cv2500_colon_InNormal <- Trans.Rang1.cv2500[c(1:39),]
Trans.Rang1.cv2500_colon_InNormal[1:9,]<- base::apply(Trans.Rang1.cv2500_colon_InNormal[1:9,], 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.cv2500_colon_InNormal[10:39,]<- base::apply(Trans.Rang1.cv2500_colon_InNormal[10:39,], 2, function(x){
  (x-min(x))/(max(x)-min(x))
})


#### 6.1.TumorOnly_Method_evaluation_Prepare.R
### 1.Loading data
Cor.merged.V8.CV2500.dataset <- readRDS("Cor.merged.V8.CV2500.dataset.rds")
Cor.merged.V8.CV2500 <- Cor.merged.V8.CV2500.dataset$Cor.merged.V8.CV2500
COAD.pheno <- Cor.merged.V8.CV2500.dataset$COAD.pheno
COAD.pheno$sampleTypes
### 2.Extract tumor only samples 
COAD.pheno.tumor <- COAD.pheno[COAD.pheno$sampleTypes == "Tumor",]
Cor.merged.V8.CV2500.tumor <- Cor.merged.V8.CV2500[,rownames(COAD.pheno.tumor)]

Cor.merged.V8.CV2500.tumor.dataset <- list(Cor.merged.V8.CV2500.tumor = Cor.merged.V8.CV2500.tumor,
                                           COAD.pheno.tumor = COAD.pheno.tumor)

saveRDS(Cor.merged.V8.CV2500.tumor, file = "Tumor.Cor.merged.V8.CV2500.dataset.rds")
### 3.Twice Normalization
## 1) First cross sample scale
#normal the correlation across samples
Trans.Rang1.CV2500.tumor<- base::apply(Cor.merged.V8.CV2500.tumor, 1, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.CV2500.tumor <- t(Trans.Rang1.CV2500.tumor)

## 2)The second Normalization
# Normalize the correlation across each cell type cluster; It's the weight of each cell type cluster
cellTypesName <- rownames(Trans.Rang1.CV2500.tumor)
#grep("Adult", cellTypesName)
Trans.Rang1.CV2500.Normalized.tumor <- Trans.Rang1.CV2500.tumor
Trans.Rang1.CV2500.Normalized.tumor[grep("Adult", cellTypesName),]<- base::apply(Trans.Rang1.CV2500.Normalized.tumor[grep("Adult", cellTypesName),], 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.Rang1.CV2500.Normalized.tumor[grep("Fetal", cellTypesName),]<- base::apply(Trans.Rang1.CV2500.Normalized.tumor[grep("Fetal", cellTypesName),], 2, function(x){
  (x-min(x))/(max(x)-min(x))
})

pheatmap::pheatmap(Trans.Rang1.CV2500.Normalized.tumor,cluster_rows = F)

## 3) save result
Trans.Rang1.CV2500.Normalized.tumor.dataset <- list(Trans.Rang1.CV2500.Normalized.tumor = Trans.Rang1.CV2500.Normalized.tumor,
                                                    COAD.pheno.tumor = COAD.pheno)
saveRDS(Trans.Rang1.CV2500.Normalized.tumor.dataset, file = "Tumor.Cor.merged.V8.CV2500.scaled.Twice.dataset.rds")




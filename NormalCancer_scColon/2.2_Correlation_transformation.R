#### Raw Correlation value transformed into correlation percent for each cell
rawCorrelation <- Cor.Res.CV8000$Cor.merged
rawCorrelation[1:5,1:5]
Trans.Rang1<- base::apply(rawCorrelation, 2, function(x){
  (x-min(x))/(max(x)-min(x))
})
Trans.sum1<- base::apply(rawCorrelation, 2, function(x){
  x/sum(x)
})

#### SC3 cluster of Cor.Res
#### Step1 prepare
library(SingleCellExperiment)
library(SC3)
library(scater)
ann <- Cor.Res.CV2000$Pheno.merged
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(Trans.Rang1),
    logcounts = as.matrix(Trans.Rang1)
  ), 
  colData = ann
)
# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
plotPCA(sce, colour_by = "Cell_info")
#### Step2 Run SC3
sce <- sc3(sce, ks = 2:4, biology = F,gene_filter = FALSE, n_cores = 15)
### colData
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
#sc3_interactive(sce)
sc3_export_results_xls(sce)
# plot PCA
plotPCA(sce, 
        colour_by = "sc3_4_clusters"
)





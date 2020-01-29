#### SC3 cluster of Cor.Res
#### Step1 prepare
library(SingleCellExperiment)
library(SC3)
library(scater)

Cor.merged.log2 <- Cor.Res.CV2000$Cor.merged
ann <- Cor.Res.CV2000$Pheno.merged
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(Cor.merged),
    logcounts = as.matrix(Cor.merged)
  ), 
  colData = ann
)
# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
plotPCA(sce, colour_by = "Cell_info")
#### Step2 Run SC3
sce <- sc3(sce, ks = 2:4, biology = F,gene_filter = FALSE, n_cores = 10)
### colData
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
#sc3_interactive(sce)
sc3_export_results_xls(sce)
# plot PCA
plotPCA(sce, 
  colour_by = "sc3_4_clusters"
)
#### Step4 Plot Functions
### Consensus Matrix
sc3_plot_consensus(sce, k = 3)
sc3_plot_consensus(
  sce, k = 3, 
  show_pdata = c(
    "Cell_info",
    "sc3_3_clusters"
  )
)
# Expression Matrix
sc3_plot_expression(sce, k = 3, show_pdatac(
  "Cell_info",
))

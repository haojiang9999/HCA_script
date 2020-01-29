BiocManager::install("scater")
#### 3.1 SC3 Input
library(SingleCellExperiment)
library(SC3)
library(scater)
head(ann)
yan[1:3, 1:3]
# create a SingleCellExperiment object
# normalised and log-transformed expression matrix, is used in the main clustering algorithm
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(yan),
    logcounts = log2(as.matrix(yan) + 1)
  ), 
  colData = ann
)
# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
isSpike(sce)
plotPCA(sce, colour_by = "cell_type1")

#### 3.2 Run SC3
sce <- sc3(sce, ks = 2:4, biology = TRUE)

sc3_interactive(sce)
sc3_export_results_xls(sce)

#### 3.3colData
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

plotPCA(
  sce, 
  colour_by = "sc3_3_clusters", 
  size_by = "sc3_3_log2_outlier_score"
)



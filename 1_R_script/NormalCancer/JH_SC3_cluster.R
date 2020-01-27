###### Function of cluster by SC3
JH_SC3_cluster <- function(expr, ann, ks=2:4,gene_filter = FALSE, n_cores = 15){
  #ann <- Pheno.merged
  #expr <- Min_max_norm.1500
  #expr <- expr[,1:100]
  #### In order to use SC3 have to remove cell had zero variance 
  ZeroVar<- which(apply(expr, 2, var)==0)
  ann <- ann[!rownames(ann) %in% names(ZeroVar), ]
  expr<- expr[,!colnames(expr) %in% names(ZeroVar)]
  sampleNames <- colnames(expr)
  #### Step1 construct sc3 object
  if (!require(SingleCellExperiment)) BiocManager::install(SingleCellExperiment)
  if (!require(SC3)) BiocManager::install(SC3)
  if (!require(scater)) BiocManager::install(scater)
  require(SingleCellExperiment)
  require(SC3)
  require(scater)
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(expr),
      logcounts = as.matrix(expr)
    ), 
    colData = ann
  )
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  #### Step2 Run SC3
  sce <- sc3(sce, ks = ks, biology = F,gene_filter = gene_filter, n_cores = n_cores)
  SC3_cluster <- metadata(sce)$sc3$consensus
  #SC3_cluster$`2`$hc$labels
  #class(SC3_cluster)
  #### Step3 Correct the labels of each cluster result
  SC3_cluster<- lapply(SC3_cluster, function(x){
    x$hc$labels <- sampleNames
    return(x)
  })
  return(SC3_cluster)
}

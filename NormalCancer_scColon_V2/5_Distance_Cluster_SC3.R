#### Distance Cluster using SC3
#### SC3 cluster of Cor.Res
#### Step1 prepare
library(SingleCellExperiment)
library(SC3)
library(scater)
ann <- Pheno.merged
expr <- Min_max_norm.1500
#### In order to use SC3 have to remove cell had zero variance 
ZeroVar<- which(apply(expr, 2, var)==0)
g[!rownames(g) %in% remove, ]
ann <- ann[!rownames(ann) %in% names(ZeroVar), ]
expr<- expr[,!colnames(expr) %in% names(ZeroVar)]
# construct sc3 object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(expr),
    logcounts = as.matrix(expr)
  ), 
  colData = ann
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
plotPCA(sce, colour_by = "ajcc_pathologic_tumor_stage")
#### Step2 Run SC3
sce <- sc3(sce, ks = 2:4, biology = F,gene_filter = FALSE, n_cores = 15)
##(optional) sc3_estimate_k
sce <- sc3_estimate_k(sce)
str(metadata(sce)$sc3)
sc3_plot_silhouette(sce, k = 4)

colData(sce)
dev.off() 
sc3_plot_consensus(
  sce, k = 3)

sc3_plot_consensus(
  sce, k = 3, 
  show_pdata = c(
    "Cell_info", 
    "CellType",
    "sc3_3_clusters"
  )
)


sc3_plot_expression(sce, k = 3,
                    show_pdata = c(
                      "Cell_info", 
                      "CellType",
                      "sc3_3_clusters"
                    ))


names(metadata(sce)$sc3$consensus)
names(metadata(sce)$sc3$consensus$`3`)
K3_cluster <- metadata(sce)$sc3$consensus$`3`

K3_cluster$hc
reorder

source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Min_max_norm.2000,show_colnames = F,
           annotation_col = Pheno.merged[,c(1,2)], cluster_cols = K3_cluster$hc)

### It looks great


#### Build T cells in colon cancer Zemin_GSE108989
### 1.Read data set
CRC.TCell.S11138.TPM.dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Zemin_GSE108989_T_cell_in_Colon_Cancer/GSE108989_CRC.TCell.S11138.TPM.dataset.rds")
CRC_TCell.S11138.TPM.exp <- CRC.TCell.S11138.TPM.dataset$GSE108989_CRC_TCell.S11138.TPM.exp
Cell_annotation <- CRC.TCell.S11138.TPM.dataset$Cell_annotation
### 2.Filter cell with low TPM sum (Remove low quality cells)
# I consider all cell pass the QC, all above 8000000
summary(colSums(CRC_TCell.S11138.TPM.exp))
Cell_annotation$cellInfo
table(Cell_annotation$cell_Info)
### 3.Build reference cell expression table
# construct _reference_panel_by_cell types
CRC.Tcell.mergeTable <- list()
CRC.Tcell.cellTypes<- as.character(unique(Cell_annotation$cell_Info))
table(Cell_annotation$cell_Info)
length()
# i="CD4+CD25hi cells from peripheral blood"
for (i in CRC.Tcell.cellTypes) {
  # select samples for each cluster
  clusterSample <- Cell_annotation[Cell_annotation$cell_Info == i, ]$UniqueCell_ID
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- CRC_TCell.S11138.TPM.exp[, clusterSample]
  #normalize expression by cell
  CRC.Tcell.mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(CRC.Tcell.mergeTable[[i]]) <- paste0("CRC_T_Cell_",i)
}

CRC.Tcell.ref <- do.call("cbind",CRC.Tcell.mergeTable)
head(CRC.Tcell.ref)
dim(CRC.Tcell.ref)
summary(colSums(CRC.Tcell.ref))
### 4.Build reference cell expression table by cluster result
## Cluster results from the paper
# construct _reference_panel_by cluster resaults
table(Cell_annotation$Cluster_anno)
CRC.Tcell.mergeTable.Cluster <- list()
CRC.Tcell.cellCluster<- as.character(unique(Cell_annotation$Cluster_anno))
# i="CD4+T_N"
for (i in CRC.Tcell.cellCluster[1:20]) {
  # select samples for each cluster
  clusterSample <- Cell_annotation[Cell_annotation$Cluster_anno %in% i, ]$UniqueCell_ID
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- CRC_TCell.S11138.TPM.exp[, clusterSample]
  #normalize expression by cell
  CRC.Tcell.mergeTable.Cluster[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(CRC.Tcell.mergeTable.Cluster[[i]]) <- paste0("CRC_T_Cell_",i)
}
CRC.Tcell.Cluster.ref <- do.call("cbind",CRC.Tcell.mergeTable.Cluster)
head(CRC.Tcell.Cluster.ref)
dim(CRC.Tcell.Cluster.ref)
summary(colSums(CRC.Tcell.Cluster.ref))


### 4.Save the data
saveRDS(CRC.Tcell.ref, file = "Zemin.CRC.Tcell.ref.rds")
saveRDS(CRC.Tcell.Cluster.ref, file = "Zemin.CRC.Tcell.Cluster.ref.rds")

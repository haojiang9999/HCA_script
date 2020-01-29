# construct _reference_panel_by_cell types
mergeTable <- list()
i=1
for (i in 1:max(as.numeric(cell.lineage$cellTypeNum))) {
  # select samples for each cluster
  clusterSample <- cell.lineage[cell.lineage$cellTypeNum == i, ]$cellName
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- FPKM[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("Embryo_",i)
}
Embryo.ref.panel.RPKM <- do.call("cbind",mergeTable)
head(Embryo.ref.panel.RPKM)
# change column names to cell type
colnames(Embryo.ref.panel.RPKM) <- paste0(colnames(Embryo.ref.panel.RPKM),"_",cellType$cellType)
head(Embryo.ref.panel.RPKM)
summary(colSums(RPKM))
#rm(Embryo.ref.panel)
## save the resault
saveRDS(Embryo.ref.panel.RPKM, file = "referPanel_of_embryo_E_MTAB_3929_RPKM.rds")

#### build Tang adult colon reference panel
#Expression matrix
filePath <- "/stor/jianghao/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE103154_adult_tissues/R_GSE103154_Tang_adult_colon_normal"
file <- paste0(filePath,"/GSE103254.adult.tpm.cellFiltered.rds")
GSE103254.adult.tpm <- readRDS(file)
exprMat <- GSE103254.adult.tpm$GSE103254.adult.cellFlitered.tpm
cellInfo <- GSE103254.adult.tpm$cellInfo.cellFiltered
hcResault.final # Author cluster cell types
hcResault.final$hc_10 <- as.numeric(hcResault.final$hc_10)
### filter reference cells sum of TPM >8000 
summary(colSums(exprMat))
keepSample <- colnames(exprMat)[colSums(exprMat) >= 800000]
mergeTable <- list()
for (i in 1:max(hcResault.final$hc_10)) {
  # select samples for each cluster
  clusterSample <- hcResault.final[hcResault.final$hc_10 == i, ]$Sample 
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- exprMat[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- unique(hcResault.final[hcResault.final$hc_10 == i, ]$cellType)
}

Tang.GI.Adult.ref.panel <- do.call("cbind",mergeTable)
head(Tang.GI.Adult.ref.panel)
colSums(Tang.GI.Adult.ref.panel)

#### Select genes in the feature sets ###
table(rowSums(Tang.GI.Adult.ref.panel) > 3)
# genes expression >3 in all cell types
Tang.GI.Adult.ref.panel.filtered <- Tang.GI.Adult.ref.panel[apply(Tang.GI.Adult.ref.panel, 1, sum) > 3, ] 
# gene expressed in more than 3 cell types
# apply(refPanel, 1, function(x) sum(x !=0) )

# Output the data
saveRDS(Tang.GI.Adult.ref.panel, file = "Tang.GI.Adult.ref.panel.rds")
saveRDS(refPanel.filtered, file = "refPanel.filtered_of_fetal.rds")
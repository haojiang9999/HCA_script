#### Build Tang Adult colon normal cell reference panel
### 1.Read data set
Tang_Adult_colon_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Tang_Colon_GI_development_dataset/Tang_Adult_colon_dataset.rds")
Adult.colon.exp.TPM <- Tang_Adult_colon_dataset$Adult.colon.exp.TPM
Adult.hcResault.paper <- Tang_Adult_colon_dataset$Adult.hcResault.paper
Adult.Gene.Anno <- Tang_Adult_colon_dataset$Adult.Gene.Anno
### 2.Filter cell with low TPM sum (Remove low quality cells)
# TPM > = 80,0000
summary(colSums(Adult.colon.exp.TPM))
Adult.keep <- colnames(Adult.colon.exp.TPM)[colSums(Adult.colon.exp.TPM) >= 800000]
# 1891 -> 1387
### 3.Build reference cell expression table
Adult.mergeTable <- list()
for (i in unique(Adult.hcResault.paper$cellType)) {
  # i = "Goblet_1"
  print(i)
  # select samples for each cluster
  clusterSample <- Adult.hcResault.paper[Adult.hcResault.paper$cellType == i, ]$Sample 
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% Adult.keep]
  #class(clusterSample)
  exprTable <- Adult.colon.exp.TPM[, clusterSample]
  #normalize expression by cell
  Adult.mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(Adult.mergeTable[[i]]) <- paste0("Adult_",i)
}
Tang.Adult.colon.ref <- do.call("cbind",Adult.mergeTable)
### Without filter genes

### 4.Save the data
saveRDS(Tang.Adult.colon.ref, file = "Tang.Adult.colon.ref.rds")

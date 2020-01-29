#### Build Tang Fetal GI normal cell reference panel
### 1.Read data set
Tang_Fetal_GI_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Tang_Colon_GI_development_dataset/Tang_Fetal_GI_dataset.rds")
Fetal.GI.exp.TPM <- Tang_Fetal_GI_dataset$Fetal.GI.exp.TPM
Fetal.Cluster.paper <- Tang_Fetal_GI_dataset$Fetal.Cluster.paper
Esophagus_KNN <- Fetal.Cluster.paper$Esophagus_KNN
LIntes_KNN <- Fetal.Cluster.paper$LIntes_KNN
SIntes_KNN <- Fetal.Cluster.paper$SIntes_KNN
Stomach_KNN <- Fetal.Cluster.paper$Stomach_KNN
Feathers_organs <- Fetal.Cluster.paper$Feathers_organs
### 2.Filter cell with low TPM sum (Remove low quality cells)
# TPM > = 80,0000
summary(colSums(Fetal.GI.exp.TPM))
Fetal.GI.exp.TPM[,colSums(Fetal.GI.exp.TPM) == 0] # Expression were all 0
Fetal.keep <- colnames(Fetal.GI.exp.TPM)[colSums(Fetal.GI.exp.TPM) >= 800000]
# 5290 -> 4748
### 3.Build reference cell expression table
###################### Esophagus ##########################
#subset cell expression data  by organs
mergeTable <- list()
for (i in 1:max(Esophagus_KNN$Cluster)) {
  # select samples for each cluster
  # i = 1
  clusterSample <- Esophagus_KNN[Esophagus_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% Fetal.keep]
  #class(clusterSample)
  exprTable <- Fetal.GI.exp.TPM[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("Esophagus_",i)
}
Fetal.Esophagus.ref <- do.call("cbind",mergeTable)
#Esophagus_KNN[,clusterSample]
#Esophagus_KNN[Esophagus_KNN$Sample %in% clusterSample,]
###################### Stomach ##########################
#subset cell expression data  by organs
mergeTable <- list()
for (i in 1:max(Stomach_KNN$Cluster)) {
  # select samples for each cluster
  # i = 1
  clusterSample <- Stomach_KNN[Stomach_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% Fetal.keep]
  #class(clusterSample)
  exprTable <- Fetal.GI.exp.TPM[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("Stomach_",i)
}
Fetal.Stomach.ref <- do.call("cbind",mergeTable)
###################### SIntes ##########################
#subset cell expression data  by organs
mergeTable <- list()
for (i in 1:max(SIntes_KNN$Cluster)) {
  # select samples for each cluster
  # i = 1
  clusterSample <- SIntes_KNN[SIntes_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% Fetal.keep]
  #class(clusterSample)
  exprTable <- Fetal.GI.exp.TPM[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("SIntes_",i)
}
Fetal.SIntes.ref <- do.call("cbind",mergeTable)
###################### LIntes ##########################
#subset cell expression data  by organs
mergeTable <- list()
for (i in 1:max(LIntes_KNN$Cluster)) {
  # select samples for each cluster
  # i = 1
  clusterSample <- LIntes_KNN[LIntes_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% Fetal.keep]
  #class(clusterSample)
  exprTable <- Fetal.GI.exp.TPM[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("LIntes_",i)
}
Fetal.LIntes.ref <- do.call("cbind",mergeTable)


#### 4.Merge all the columns together
Fetal.GI.ref <- cbind(Fetal.Esophagus.ref, Fetal.Stomach.ref, Fetal.SIntes.ref, Fetal.LIntes.ref)
colSums(Fetal.GI.ref)
# change column names to cell type
colnames(Fetal.GI.ref) <- paste0(colnames(Fetal.GI.ref),"_",Feathers_organs$Cell.Type)
View(Fetal.GI.ref)
### Without filter genes
### 5.Save the data
saveRDS(Fetal.GI.ref, file = "Tang.Fetal.GI.ref.rds")


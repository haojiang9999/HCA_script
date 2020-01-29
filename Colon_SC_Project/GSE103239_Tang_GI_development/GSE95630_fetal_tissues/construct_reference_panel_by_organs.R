# reference pannel construction 
###################### Esophagus ##########################
#subset cell expression data  by organs
Esophagus <- fetal.tissues[,as.character(Esophagus_KNN$Sample)] # Caution：convert factor to character
summary(colSums(Esophagus)) # NOT all the UMIs was known genes
# sample filtered by total TPM 800000
keepSample <- colnames(Esophagus)[colSums(Esophagus) >= 800000]
mergeTable <- list()
for (i in 1:max(Esophagus_KNN$Cluster)) {
  # select samples for each cluster
  clusterSample <- Esophagus_KNN[Esophagus_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- Esophagus[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("Esophagus_",i)
}
Esophagus.ref.panel <- do.call("cbind",mergeTable)
head(Esophagus.ref.panel)

###################### Stomach ##########################
#subset cell expression data  by organs
Stomach <- fetal.tissues[,as.character(Stomach_KNN$Sample)] # Caution：convert factor to character
summary(colSums(Stomach)) # NOT all the UMIs was known genes
# sample filtered by total TPM 800000
keepSample <- colnames(Stomach)[colSums(Stomach) >= 800000]
mergeTable <- list()
for (i in 1:max(Stomach_KNN$Cluster)) {
  # select samples for each cluster
  clusterSample <- Stomach_KNN[Stomach_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- Stomach[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("Stomach_",i)
}
Stomach.ref.panel <- do.call("cbind",mergeTable)
head(Stomach.ref.panel)

###################### SIntes ##########################
#subset cell expression data  by organs
SIntes <- fetal.tissues[,as.character(SIntes_KNN$Sample)] # Caution：convert factor to character
summary(colSums(SIntes)) # NOT all the UMIs was known genes
# sample filtered by total TPM 800000
keepSample <- colnames(SIntes)[colSums(SIntes) >= 800000]
mergeTable <- list()
for (i in 1:max(SIntes_KNN$Cluster)) {
  # select samples for each cluster
  clusterSample <- SIntes_KNN[SIntes_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- SIntes[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("SIntes_",i)
}
SIntes.ref.panel <- do.call("cbind",mergeTable)
head(SIntes.ref.panel)

###################### LIntes ##########################
#subset cell expression data  by organs
LIntes <- fetal.tissues[,as.character(LIntes_KNN$Sample)] # Caution：convert factor to character
summary(colSums(LIntes)) # NOT all the UMIs was known genes
# sample filtered by total TPM 800000
keepSample <- colnames(LIntes)[colSums(LIntes) >= 800000]
mergeTable <- list()
for (i in 1:max(LIntes_KNN$Cluster)) {
  # select samples for each cluster
  clusterSample <- LIntes_KNN[LIntes_KNN$Cluster == i, ]$Sample
  clusterSample <- as.character(clusterSample)
  clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- LIntes[, clusterSample]
  #normalize expression by cell
  mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(mergeTable[[i]]) <- paste0("LIntes_",i)
}
LIntes.ref.panel <- do.call("cbind",mergeTable)
head(LIntes.ref.panel)
# merge all the columns together
refPanel <- cbind(Esophagus.ref.panel, Stomach.ref.panel, SIntes.ref.panel, LIntes.ref.panel)
colSums(refPanel)
# change column names to cell type
colnames(refPanel) <- paste0(colnames(refPanel),"_",Feathers_organs$Cell.Type)
head(refPanel)
colnames(refPanel)
summary(colSums(refPanel))

#### Select genes in the feature sets ###
table(rowSums(refPanel) > 3)
# genes expression >3 in all cell types
refPanel.filtered <- refPanel[apply(refPanel, 1, sum) > 3, ] 
# gene expressed in more than 3 cell types
# apply(refPanel, 1, function(x) sum(x !=0) )

# Output the data
saveRDS(refPanel, file = "referPanel_of_fetal.rds")
saveRDS(refPanel.filtered, file = "refPanel.filtered_of_fetal.rds")

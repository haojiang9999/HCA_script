#### Build_RCA_cell_scReference.R
### 1.Read data set
RCA_all_COUNT_normalized_uniGene_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_all_COUNT_normalized_uniGene_dataset.rds")
NM.all.COUNT.normalized.uniGene <- RCA_all_COUNT_normalized_uniGene_dataset$NM.all.COUNT.normalized.uniGene
Tumor.all.COUNT.normalized.uniGene <- RCA_all_COUNT_normalized_uniGene_dataset$Tumor.all.COUNT.normalized.uniGene
NM.all.phenoType <- RCA_all_COUNT_normalized_uniGene_dataset$NM.all.phenoType
tumor.all.phenoType <- RCA_all_COUNT_normalized_uniGene_dataset$tumor.all.phenoType
### 2.Filter cell with low TPM sum (Remove low quality cells)
# I consider all cell pass the QC, all above 8000000
summary(colSums(NM.all.COUNT.normalized.uniGene))
summary(colSums(Tumor.all.COUNT.normalized.uniGene))
table(NM.all.phenoType$All_cellTypes)
table(tumor.all.phenoType$All_cellTypes)
### 3.Build reference cell expression table
# construct _reference_panel_by_cell types
RCA.nonEpi.mergeTable.NM <- list()
table(NM.all.phenoType$All_cellTypes)
NM.all.phenoType$NM_ID
table(tumor.all.phenoType$All_cellTypes)
# i="Bcell"
############ Normal #############
for (i in c("Bcell", "Endothelial", "Fibroblast","Macrophage","MastCell","Tcell")) {
  # select samples for each cluster
  clusterSample <- rownames(NM.all.phenoType[NM.all.phenoType$All_cellTypes == i, ])
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- NM.all.COUNT.normalized.uniGene[, clusterSample]
  #normalize expression by cell
  RCA.nonEpi.mergeTable.NM[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(RCA.nonEpi.mergeTable.NM[[i]]) <- paste0("RCA_Normal_",i)
}
RCA.nonEpi.ref.NM <- do.call("cbind",RCA.nonEpi.mergeTable.NM)
head(RCA.nonEpi.ref.NM)

############ Tumor #############
RCA.nonEpi.mergeTable.Tumor <- list()
#### Because tumor had only one MastCell So I remove it

for (i in c("Bcell", "Endothelial", "Fibroblast","Macrophage","Tcell")) {
  # select samples for each cluster
  clusterSample <- rownames(tumor.all.phenoType[tumor.all.phenoType$All_cellTypes == i, ])
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- Tumor.all.COUNT.normalized.uniGene[, clusterSample]
  #normalize expression by cell
  RCA.nonEpi.mergeTable.Tumor[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(RCA.nonEpi.mergeTable.Tumor[[i]]) <- paste0("RCA_Tumor_",i)
}
RCA.nonEpi.ref.Tumor <- do.call("cbind",RCA.nonEpi.mergeTable.Tumor)
head(RCA.nonEpi.ref.Tumor)

### 4.Save data
RCA.nonEpi.ref <- cbind(RCA.nonEpi.ref.NM,RCA.nonEpi.ref.Tumor)
head(RCA.nonEpi.ref)
saveRDS(RCA.nonEpi.ref, file = "RCA.CRC.nonEpi.ref.rds")






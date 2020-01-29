#### Build Lanner Pre-implantation Embryo reference panel
### 1.Read data set
Pre_implantation_E_MTAB_3929_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Preimplantation_Embryos_Lanner/Pre_implantation_E_MTAB_3929_dataset.rds")
Preim.embryo.exp.COUNT.normal <- Pre_implantation_E_MTAB_3929_dataset$Preim.embryo.exp.COUNT.normal
Preim.embryo.pheno <- Pre_implantation_E_MTAB_3929_dataset$Preim.embryo.pheno
### 2.Filter cell with low TPM sum (Remove low quality cells)
# I consider all cell pass the QC
### 3.Build reference cell expression table
# construct _reference_panel_by_cell types
Preim.embryo.mergeTable <- list()
table(Preim.embryo.pheno$CellTypes)
# i="Prelineage"
Preim.embryo.pheno$Source.Name
for (i in c("Prelineage", "epiblast", "primitive endoderm")) {
  # select samples for each cluster
  clusterSample <- Preim.embryo.pheno[Preim.embryo.pheno$CellTypes == i, ]$Source.Name
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- Preim.embryo.exp.COUNT.normal[, clusterSample]
  #normalize expression by cell
  Preim.embryo.mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(Preim.embryo.mergeTable[[i]]) <- paste0("Pre_Embryo_",i)
}
## Seperate trophectoderm cells into mural and polar cell types
for (i in c("mural", "polar")) {
  # select samples for each cluster
  clusterSample <- Preim.embryo.pheno[Preim.embryo.pheno$Characteristics.inferred.trophectoderm.subpopulation. == i, ]$Source.Name
  clusterSample <- as.character(clusterSample)
  #clusterSample <- clusterSample[clusterSample %in% keepSample]
  #class(clusterSample)
  exprTable <- Preim.embryo.exp.COUNT.normal[, clusterSample]
  #normalize expression by cell
  Preim.embryo.mergeTable[[i]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
  names(Preim.embryo.mergeTable[[i]]) <- paste0("Pre_Embryo_trophectoderm_",i)
}
Preim.embryo.ref <- do.call("cbind",Preim.embryo.mergeTable)
head(Preim.embryo.ref)

### 4.Save the data
saveRDS(Preim.embryo.ref, file = "Lanner.Preim.embryo.ref.rds")


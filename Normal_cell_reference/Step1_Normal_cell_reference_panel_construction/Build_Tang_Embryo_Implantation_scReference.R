#### Build Tang Embryo implantation reference panel
### 1.Read data set
Tang_Normal_embryos_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Tang_human_implantation_GSE109555/Tang_Normal_embryos_dataset.rds")
Tang.normal.embryo.exp.TPM <- Tang_Normal_embryos_dataset$Tang.normal.embryo.exp.TPM
Tang.normal.embryo.Sample_Information <- Tang_Normal_embryos_dataset$Tang.normal.embryo.Sample_Information
Tang.embryo.Gene.Anno <- Tang_Normal_embryos_dataset$Tang.embryo.Gene.Anno
### 2.Filter cell with low TPM sum (Remove low quality cells)
# TPM > = 80,0000
summary(colSums(Tang.normal.embryo.exp.TPM))
# The lowest was 99,0000 so no one was removed
### 3.Build reference cell expression table
table(Tang.normal.embryo.Sample_Information$Day)
table(Tang.normal.embryo.Sample_Information$Lineage)
Tang.normal.embryo.mergeTable <- list()
for (Day in c("D6","D8","D10","D12")) {
  # select samples for each cluster
  # Day = "D8"
  DaySample <- Tang.normal.embryo.Sample_Information[Tang.normal.embryo.Sample_Information$Day == Day, ]
  for (lineage in c("EPI", "PE", "TE")){
    # lineage = "TE"
    clusterSample <- DaySample[DaySample$Lineage == lineage, ]$Sample
    clusterSample <- as.character(clusterSample)
    #clusterSample <- clusterSample[clusterSample %in% Fetal.keep]
    exprTable <- Tang.normal.embryo.exp.TPM[, clusterSample]
    #normalize expression by cell
    CellName <- paste0("Tang_embryo_",Day,"_",lineage)
    Tang.normal.embryo.mergeTable[[CellName]] <- as.data.frame(rowSums(exprTable)/length(clusterSample))
    names(Tang.normal.embryo.mergeTable[[CellName]]) <- CellName
  }
}
Tang.normal.embryo.ref <- do.call("cbind",Tang.normal.embryo.mergeTable)
#### Test some marker genes
# High in EPI
GeneName <- c("DPPA5", "KHDC3L", "TDGF1", "NANOG", "SOX2", "CXCL12", "PRDM14", "THY1")
Tang.normal.embryo.ref[GeneName,]
# Looks it fits the paper resaults
### 4.Save the data
saveRDS(Tang.normal.embryo.ref, file = "Tang.Normal.embryo.ref.rds")









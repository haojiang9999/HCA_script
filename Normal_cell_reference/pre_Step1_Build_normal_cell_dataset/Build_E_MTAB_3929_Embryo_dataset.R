#### Build_E_MTAB_3929_Embryo_dataset ####
filePath.counts <- "/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/counts.txt"
filePath.rpkm <- "/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/rpkm.txt"
E_MTAB_3929_Embryo.counts <- read.table(filePath.counts,header = T)
E_MTAB_3929_Embryo.rpkm <- read.table(filePath.rpkm,header = T)
head(E_MTAB_3929_Embryo.counts[1:10,1:10])
summary(colSums(E_MTAB_3929_Embryo.counts)) 
summary(colSums(E_MTAB_3929_Embryo.rpkm)) 
hist(as.numeric(E_MTAB_3929_Embryo.rpkm[2,]))
hist(log2(as.numeric(E_MTAB_3929_Embryo.rpkm[2,])))
hist(as.numeric(E_MTAB_3929_Embryo.counts[2,]))
hist(log2(as.numeric(E_MTAB_3929_Embryo.counts[2,])))
# the resault showed this was no-transformed tpm value
# But some samples sum was too low not filtered

### Step2 Phenotype information
## extract from supplementary data
pd <- read.delim2("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/E-MTAB-3929.sdrf.txt")
## loading cell label from chen wenchang
load("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/em_label1.RData")
load("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/em_label2.RData")
label2$label.lng[1,]
cell.lineage <- data.frame(cellName = colnames(E_MTAB_3929_Embryo.rpkm), 
                           cellTypeNum = as.factor(label2$label.lng[1,]))
cellType <- data.frame(cellTypeNum = c("1","2","3","4"), 
                      cellType = c("EPI","Prelineage","PE","TE"))
cell.lineage <- dplyr::left_join(cell.lineage, cellType, by = "cellTypeNum")

### Step3 Gene features
# they don`t provide

### Step4 Build datasets
E_MTAB_3929_Embryo_datasets <- list(E_MTAB_3929_Embryo.counts = E_MTAB_3929_Embryo.counts,
                                    E_MTAB_3929_Embryo.rpkm = E_MTAB_3929_Embryo.rpkm,
                                            pd = pd,
                                    cell.lineage = cell.lineage)
saveRDS(E_MTAB_3929_Embryo_datasets, file = "E_MTAB_3929_Embryo_datasets.rds")







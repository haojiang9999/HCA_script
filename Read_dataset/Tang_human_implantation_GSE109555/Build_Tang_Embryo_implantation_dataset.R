#### Build Tang Embryo implantation dataset
########################## Adult colon tissue #############################################
## 1.Read expression data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109555
# read expression table
all.TPM.file <- "/stor/jianghao/GEO/Normal_tissues/GSE109555_Tang_human_implantation/GSE109555_All_Embryo_TPM.txt.gz" 
All_Embryo_TPM <- read.table(all.TPM.file)
head(All_Embryo_TPM[1:10.1:10])
# barcod files dose not work
GSE109555.RAW.file <- "/stor/jianghao/GEO/Normal_tissues/GSE109555_Tang_human_implantation/GSE109555_RAW.tar"
GSE109555_RAW <- read.delim(GSE109555.RAW.file)
# Cells using TrioSeq tech
TrioSeq.DataInfo.file <- "/stor/jianghao/GEO/Normal_tissues/GSE109555_Tang_human_implantation/GSE109555_TrioSeq_DataInfo.txt.gz"
TrioSeq_DataInfo <- read.table(TrioSeq.DataInfo.file)
# TrioSeq cell expression
TrioSeq.TPM.file <- "/stor/jianghao/GEO/Normal_tissues/GSE109555_Tang_human_implantation/GSE109555_TrioSeq_TPM.txt.gz"
TrioSeq_TPM <- read.table(TrioSeq.TPM.file)
############# All_Embryo_TPM not using TrioSeq tech and using paper used normal ones ###################
## 2.Original paper cluster and celltype annotation
### read metadata of embryos 
library(readr)
Embryo_statistics <- read_csv("Embryo_statistics.csv")
View(Embryo_statistics)
## Find embryos had three major lineages which were considered normal embryos
Normal.embryos <- Embryo_statistics[Embryo_statistics$Normal == "Yes", ]$Sample
### Sample information from paper
# Reconstituting the transcriptome and DNA methylome landscapes of human implantation
Sample_Information <- read_csv("Supplementary Table 2 Sample Information.csv")
table(Sample_Information$Ori_Day_Emb %in% Normal.embryos)
Sample_Information.Normal <- Sample_Information[Sample_Information$Ori_Day_Emb %in% Normal.embryos,]
Samples.Normal <- Sample_Information.Normal$Sample
Tang.normal.embryo.TPM <- All_Embryo_TPM[,Samples.Normal]
summary(colSums(Tang.normal.embryo.TPM))
### 3.Gene annotation
Tang.embryo.geneNames <- rownames(Tang.normal.embryo.TPM)
source("/data8t_4/JH/MyJobs/1_R_script/GeneSymbol2GeneID.R")
Tang.embryo.Gene.Anno <- GeneSymbol2GeneID(as.character(Tang.embryo.geneNames), toType = c("ENSEMBL","ENTREZID"))
length(unique(Tang.embryo.Gene.Anno$SYMBOL))
## 4.Build Adult colon dataset
Tang_Normal_embryos_dataset <- list(Tang.normal.embryo.exp.TPM = Tang.normal.embryo.TPM, Tang.normal.embryo.Sample_Information = Sample_Information.Normal ,
                                    Tang.normal.embryo.statistics = Embryo_statistics, Tang.embryo.Gene.Anno = Tang.embryo.Gene.Anno)


#### Save the data set ####
saveRDS(Tang_Normal_embryos_dataset, file = "Tang_Normal_embryos_dataset.rds")




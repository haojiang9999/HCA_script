### read all data set
RCA_all_FPKM_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_all_FPKM_dataset.rds")
NM.all.fpkm.exp <- RCA_all_FPKM_dataset$NM.all.fpkm.exp
tumor.all.fpkm.exp <- RCA_all_FPKM_dataset$tumor.all.fpkm.exp
NM.all.phenoType <- RCA_all_FPKM_dataset$NM.all.phenoType
tumor.all.phenoType <- RCA_all_FPKM_dataset$tumor.all.phenoType
#### expression data normalization Using pQ
library(devtools)
load_all("./NODES_0.0.0.9010/NODES/")
all.fpkm.exp <- cbind(NM.all.fpkm.exp, tumor.all.fpkm.exp)
all.fpkm.exp.pQ <- pQ(data = all.fpkm.exp)
saveRDS(all.fpkm.exp.pQ, file = "RCA_all.fpkm.exp.pQ.rds")
# low expression geen was 31.58095
head(epithelial.fpkm.exp.stemTA.pQ)
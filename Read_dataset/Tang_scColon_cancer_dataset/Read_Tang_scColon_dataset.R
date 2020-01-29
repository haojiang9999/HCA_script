#### Tang scColon cancer dataset building
## Folder of processed data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/GSE97693_processed/"
fileNames <- list.files(path=filePath,
                        pattern="txt.gz$")
files <- paste0(filePath,fileNames)
# Find the colon pateint CRC files by FPKM and TPM 
filesFPKM <- fileNames[grepl("(?=.*FPKM)(?=.*CRC)", fileNames, perl = T)]
filesTPM <- fileNames[grepl("(?=.*TPM)(?=.*CRC)", fileNames, perl = T)]
# gzfile() open .gz files
geneNamesFPKM <- read.table(gzfile(paste0(filePath,filesFPKM)[1]), header=FALSE, sep="\t")[,1]     # gene names
geneNamesTPM <- read.table(gzfile(paste0(filePath,filesTPM)[1]), header=FALSE, sep="\t")[,1] 
### 1.Read expression files
#Using multiple core to speed up
library(parallel)
#find cores number
detectCores(logical = F)
#set core numbers
mc <- getOption("mc.cores", 10)
##cbind() will change the data some time so cbind.data.frame maybe help
df.FPKM<- do.call(cbind.data.frame,mclapply(paste0(filePath,filesFPKM), 
                                            function(fn)read.table(gzfile(fn),header=F, sep="\t")[,2],mc.cores = mc))
class(df.FPKM)
df.TPM<- do.call(cbind.data.frame,mclapply(paste0(filePath,filesTPM), 
                                           function(fn)read.table(gzfile(fn),header=F, sep="\t")[,2],mc.cores = mc))
#Stop
stopCluster(mc)
# rownames of gene names
rownames(df.FPKM) <- geneNamesFPKM
rownames(df.TPM) <- geneNamesTPM
# columns of sample names
colnames(df.FPKM) <- sub(".FPKM.txt.gz","",filesFPKM)
colnames(df.TPM) <- sub(".TPM.txt.gz","",filesTPM)
colnames(df.TPM) <- sub("_scTrioSeq2Rna_scTrioSeq2Rna_","_scTrioSeq2Rna_",colnames(df.TPM))
### 2.phenoType data 
# df.FPKM
GEOid<- unlist(lapply(strsplit(colnames(df.FPKM),"_"), '[[', 1))
patientID <- unlist(lapply(strsplit(colnames(df.FPKM),"_"), '[[', 3))
sampleRegion <- unlist(lapply(strsplit(colnames(df.FPKM),"_"), '[[', 4))
cellNum <- unlist(lapply(strsplit(colnames(df.FPKM),"_"), '[[', 5))
FPKM.500.pheno<- data.frame(GEOid = GEOid, patientID = patientID, sampleRegion = sampleRegion,
                            cellNum = cellNum)
rownames(FPKM.500.pheno) <- colnames(df.FPKM)
# df.TPM
#TPM.688.pheno
GEOid<- unlist(lapply(strsplit(colnames(df.TPM),"_"), '[[', 1))
patientID <- unlist(lapply(strsplit(colnames(df.TPM),"_"), '[[', 3))
sampleRegion <- unlist(lapply(strsplit(colnames(df.TPM),"_"), '[[', 4))
cellNum <- unlist(lapply(strsplit(colnames(df.TPM),"_"), '[[', 5))
TPM.688.pheno<- data.frame(GEOid = GEOid, patientID = patientID, sampleRegion = sampleRegion,
                            cellNum = cellNum)
rownames(TPM.688.pheno) <- colnames(df.TPM)
### 3.Gene annotation
source("/data8t_4/JH/MyJobs/1_R_script/GeneSymbol2GeneID.R")
library(org.Hs.eg.db)
#$BiocManager::install("AnnotationDbi")
keytypes(org.Hs.eg.db)
Gene.Anno.FPKM <- GeneSymbol2GeneID(as.character(geneNamesFPKM), toType = c("ENSEMBL","ENTREZID"))
Gene.Anno.TPM <- GeneSymbol2GeneID(as.character(geneNamesTPM), toType = c("ENSEMBL","ENTREZID"))
length(unique(Gene.Anno.FPKM$SYMBOL))
length(unique(Gene.Anno.TPM$SYMBOL))
### 4.Build dataset
Tang_scColon_dataset <- list(FPKM.500.exp = df.FPKM, FPKM.500.pheno = FPKM.500.pheno,
                             FPKM.500.Gene.Anno = Gene.Anno.FPKM, TPM.688.exp = df.TPM, 
                             TPM.688.pheno = TPM.688.pheno, TPM.688.Gene.Anno = Gene.Anno.TPM)
# save data 
saveRDS(Tang_scColon_dataset, file = "Tang_scColon_dataset.rds")






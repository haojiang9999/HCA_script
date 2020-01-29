##Merge featureCoubts output
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/GSE97693_processed/"

fileNames <- list.files(path=filePath,
                    pattern="txt.gz$")
files <- paste0(filePath,fileNames)
# seperate colon pateint CRC files by FPKM and TPM 
filesFPKM <- fileNames[grepl("(?=.*FPKM)(?=.*CRC)", fileNames, perl = T)]
filesTPM <- fileNames[grepl("(?=.*TPM)(?=.*CRC)", fileNames, perl = T)]
# gzfile() open .gz files
geneNamesFPKM <- read.table(gzfile(paste0(filePath,filesFPKM)[1]), header=FALSE, sep="\t")[,1]     # gene names
geneNamesTPM <- read.table(gzfile(paste0(filePath,filesTPM)[1]), header=FALSE, sep="\t")[,1] 
#Using multiple core to speed up
library(parallel)
#find cores number
detectCores(logical = F)
#set core numbers
mc <- getOption("mc.cores", 10)
#Using the mclapply
##cbind() will change the data some time so cbind.data.frame maybe help
df.FPKM<- do.call(cbind.data.frame,mclapply(paste0(filePath,filesFPKM), 
                                            function(fn)read.table(gzfile(fn),header=F, sep="\t")[,2],mc.cores = mc))
class(df.FPKM)
df.TPM<- do.call(cbind.data.frame,mclapply(paste0(filePath,filesTPM), 
                                           function(fn)read.table(gzfile(fn),header=F, sep="\t")[,2],mc.cores = mc))
#Stop
stopCluster(mc)
# rownames of gene neames
rownames(df.FPKM) <- geneNamesFPKM
rownames(df.TPM) <- geneNamesTPM
# columns of sample names
colnames(df.FPKM) <- sub(".FPKM.txt.gz","",filesFPKM)
colnames(df.TPM) <- sub(".FPKM.txt.gz","",filesTPM)
# save data 
saveRDS(df.FPKM, file = "GSE97693_Tang_FPKM_cells.rds")
saveRDS(df.TPM, file = "GSE97693_Tang_TPM_cells.rds")

# separate normal and tumor samples
colnames(df.FPKM)
grep("NC", colnames(df.FPKM))
colnames(df.FPKM)[grep("NC", colnames(df.FPKM))]
table(grepl("NC", colnames(df.TPM)))
table(grepl("NC", colnames(df.FPKM)))









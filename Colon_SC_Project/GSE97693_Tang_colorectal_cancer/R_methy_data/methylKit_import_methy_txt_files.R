# read the CpG sites first
library(methylKit)
# sample.list
sample.listPath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/myjob/MetSample.list"
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/GSE97693_Met"
sample.list <- read.table(sample.listPath)
filesPath.all <- paste0(filePath,'/',sample.list[,1],".CpG.txt.gz")
filesPath.list.all<- as.list(filesPath.all)
length(sample.list[[1]])
# cut the sample ID names to CRC13_LN2_377 type
sample.id <- sub("^[^C]*", "", as.character(sample.list[,1]))
# group cell by pateint ID
unlist(lapply(strsplit(sample.id,"_"), '[[', 1))
table(unlist(lapply(strsplit(sample.id,"_"), '[[', 1)))
unlist(lapply(strsplit(sample.id,"_"), '[[', 1))
treatment<-c(rep(1,409),rep(2,59),rep(4,93),rep(9,27),rep(10,124),rep(11,254),
             rep(12,23),rep(13,210),rep(14,36),rep(15,60))
length(treatment)
# read all Methy data from Tang colon cancer
objDB.all=methRead(filesPath.list.all,sample.id=as.list(sample.id),
                   treatment=treatment, # group samples by patients
                   mincov = 1, # for single cells read coverage set to 1
                   assembly="hg19",header=F, context="CpG", resolution="base",
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=5,strand.col=4,freqC.col=8 ))
format(object.size(objDB.all) , units = "MB")

getMethylationStats(objDB.all)
#3.1 Merging samples
meth.all=unite(objDB.all, destrand=FALSE)
getCorrelation(meth.all,plot=TRUE)
#3.3 Clustering samples
clusterSamples(meth.all, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth.all)
dim(meth.all)

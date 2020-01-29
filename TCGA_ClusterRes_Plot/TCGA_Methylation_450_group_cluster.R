#### Methylation Beta-value
### DE gene cluster resaults
cancerType = "COAD"
COAD.sur.df  # cluster Result
#### 1.Read TCGA dataset ####
#print = TRUE
#geneName = "ENSG00000232956"
#TCGA_PANCAN_Methylation450_dataset$Methyl450_hg19_GPL16304
#"/stor/jianghao/Myjobs/R/Read_data_set/TCGA_DNA_mythelation/COAD_PANCAN_Methylation450_dataset.rds"
filePath <- paste0("/stor/jianghao/Myjobs/R/Read_data_set/TCGA_DNA_mythelation/",cancerType,"_PANCAN_Methylation450_dataset.rds")
TCGA_PANCAN_Methylation450_dataset <- readRDS(filePath)
TCGA.Methylation450 <- TCGA_PANCAN_Methylation450_dataset[[paste0(cancerType,"Methylation450")]]
TCGA.pheno <- TCGA.Toil.norm_count.dataset[[paste0(cancerType,".pheno")]]
Methyl450_hg19_GPL16304 <- TCGA_PANCAN_Methylation450_dataset$Methyl450_hg19_GPL16304
#### 2. Select tumor samples
## sample type ###
sampleName <- colnames(TCGA.Methylation450) 
samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
table(samplType)
## "01" was primary tumor, "11' was normal tissue
tumorIndex <- samplType == "01"
#normalIndex <- samplType == "11"
tumorSamples <- sampleName[tumorIndex]
TCGA.pheno.tumor <- TCGA.pheno[tumorSamples,]
TCGA.Methylation450.Tumor <- TCGA.Methylation450[,tumorSamples]
summary(rowSums(TCGA.Methylation450.Tumor))
### data clean
TCGA.Methylation450.Tumor[is.na(TCGA.Methylation450.Tumor)] <- 0
source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
topcv <- TopCV(TCGA.Methylation450.Tumor, TopN = 2000 )

rowSums(TCGA.Methylation450.Tumor)
test1 <- TCGA.Methylation450.Tumor[rowSums(TCGA.Methylation450.Tumor) > 100,]


TCGA.Methylation450.Tumor["cg00050872",]
summary(TCGA.Methylation450.Tumor[,1:10])
#summary(TCGA.RSEM.gene.norm_count.Tumor)
library(pheatmap)
#expr <- log2(expr + 0.01)
y <- COAD.sur.df.clusterRes[rownames(COAD.sur.df.clusterRes) %in% colnames(topcv),]
topcv <- topcv[,rownames(y)]
cbind(colnames(topcv), rownames(COAD.sur.df.clusterRes))
sd
topcv

pheatmap(topcv,show_rownames=F,cluster_cols = hcTumor.cv10,
         show_colnames = F,annotation_col = y)

### sd
TopSD <- function(df, TopN = 10, MARGIN = 1){
  FeatureCV <- apply(df, MARGIN , function(x) {sd(x)})
  Feature.TopN <- names(head(sort(FeatureCV, decreasing = T), TopN))
  df[Feature.TopN,]
}
TopSD
topsd <- TopSD(TCGA.Methylation450.Tumor.clean, TopN = 1500 )
topsd
pheatmap(topsd,show_rownames=F,cluster_cols = hcTumor.cv10,
         show_colnames = F,annotation_col = y)

topsd.2 <- topsd[,rownames(COAD.sur.df.clusterRes)]
colnames(topsd) %in% rownames(y)
topsd.2 <- topsd[,rownames(COAD.sur.df.clusterRes)]



pheatmap(topsd,show_rownames=F,cluster_cols = hcTumor.cv10,
         show_colnames = F,annotation_col = COAD.sur.df.clusterRes)

table(colnames(topsd) %in% rownames(COAD.sur.df.clusterRes))
errroID1 <- colnames(topsd)[!colnames(topsd) %in% rownames(COAD.sur.df.clusterRes)]
errroID2 <- rownames(COAD.sur.df.clusterRes)[!rownames(COAD.sur.df.clusterRes) %in% colnames(topsd)]
commonID <- colnames(topsd)[colnames(topsd) %in% rownames(COAD.sur.df.clusterRes)]
topsd2 <- topsd[,c(commonID,errroID1)]
COAD.sur.df.clusterRes2 <- COAD.sur.df.clusterRes[c(commonID,errroID2),]

pheatmap(topsd2,show_rownames=F,cluster_cols = hcTumor.cv10,
         show_colnames = F,annotation_col = COAD.sur.df.clusterRes2)
x<-topsd
y<-COAD.sur.df.clusterRes

colnames(x)[colnames(x) %in% errroID1] <- as.character(paste0("a",seq(1,17)))
rownames(y)[rownames(y) %in% errroID2] <- as.character(paste0("a",seq(1,17)))

x2 <- x[,rownames(y)]

pheatmap(x2,show_rownames=F,cluster_cols = hcTumor.cv10,
         show_colnames = F,annotation_col = y)

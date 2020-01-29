### DE gene cluster resaults
cancerType = "COAD"
COAD.sur.df  # cluster Result
#### 1.Read TCGA dataset ####
#print = TRUE
#geneName = "ENSG00000232956"
filePath <- paste0("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_Pan_TOIL_RSEM_norm_count/",cancerType,"_UCSC_Toil_norm_count_dataset.rds")
TCGA.Toil.norm_count.dataset <- readRDS(filePath)
TCGA.RSEM.gene.norm_count <- TCGA.Toil.norm_count.dataset[[paste0(cancerType,".RSEM.gene.norm_count_round")]]
TCGA.pheno <- TCGA.Toil.norm_count.dataset[[paste0(cancerType,".pheno")]]
#### 2. Select tumor samples
## sample type ###
sampleName <- rownames(TCGA.pheno)
samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
table(samplType)
## "01" was primary tumor, "11' was normal tissue
tumorIndex <- samplType == "01"
#normalIndex <- samplType == "11"
TCGA.pheno.tumor <- TCGA.pheno[tumorIndex,]
tumorSamples <- rownames(TCGA.pheno.tumor)
TCGA.RSEM.gene.norm_count.Tumor <- TCGA.RSEM.gene.norm_count[,tumorSamples]
summary(TCGA.RSEM.gene.norm_count.Tumor)
TCGA.RSEM.gene.norm_count.Tumor[1:5,1:5]
#### 3.Find Top differentially expressd genes between groups
# https://support.bioconductor.org/p/91054/
groupNum <- length(unique(COAD.sur.df$COAD.cluster)) # how many groups

#library(curatedCRCData)
#data(TCGA.RNASeqV2_eset)
#targets <- pData(TCGA.RNASeqV2_eset)l
library(limma)
tumor.filtered <- TCGA.RSEM.gene.norm_count.Tumor[rowSums(TCGA.RSEM.gene.norm_count.Tumor) > 10,]
trans <- log2(tumor.filtered + 0.01)
y <- normalizeQuantiles(trans)
samplesID <- colnames(y)
type <- factor(COAD.sur.df[samplesID,]$COAD.cluster)
table(type)
type
keep <- rowSums(y > log2(11)) > 20
table(keep)
y2 <- y[keep,]
design <- model.matrix(~type)
fit <- lmFit(y2,design)
fit <- eBayes(fit,robust=TRUE,trend=TRUE)
topDEgenes <- topTable(fit,coef=2,number=500)
topDEgenes <- topDEgenes[topDEgenes$AveExpr > 0,]
rownames(topDEgenes)
#plotMDS(y2, label=type)
plotSA(fit)

#### 4.Top DE genes plot
topGenes<- rownames(topDEgenes)
expr <- y2[topGenes,]
#expr <- expr - rowMeans(expr)
expr <- t(scale(t(expr), scale = F))
#expr <- log2(expr +1)
saveRDS(expr, file = "expr_test1.rds")
library(pheatmap)
pheatmap(expr,scale = "none",show_rownames=F,cluster_cols = F,
         show_colnames = F,annotation_col = COAD.sur.df)
expr <- expr[,rownames(COAD.sur.df.clusterRes)]
#expr <- log2(expr + 0.01)
pheatmap(expr,show_rownames=F,cluster_cols = hcTumor.cv10,
         show_colnames = F,annotation_col = COAD.sur.df.clusterRes)





### 5.GGplot_COAD_RNA_Final.R
COAD.RNA.final.mx <- COAD_PanCancerAtlas_Publish_RNA_Final_dataset$COAD.RNA.final.mx
COAD.RNA.final.mx[1:5,1:5]
dim(COAD.RNA.final.mx)
hist(COAD.RNA.final.mx)
boxplot(COAD.RNA.final.mx)
## 1) Remove gene names were ? and expression 
rownames(COAD.RNA.final.mx)
COAD.RNA.final.sub <- COAD.RNA.final.mx[-c(1:29),]
COAD.RNA.final.sub[1:5,1:5]
table(is.na(COAD.RNA.final.sub))
table(rowSums(is.na(COAD.RNA.final.sub)) > 0)
COAD.RNA.final.sub[rowSums(is.na(COAD.RNA.final.sub)) > 0,]
## 2) Limma DE gene analysis
library(edgeR)
# DGEList object using the edgeR package
dge <- DGEList(counts=COAD.RNA.final.sub)
## filter

keep <- filterByExpr(dge, design.expr)
dge <- dge[keep,,keep.lib.sizes=FALSE]

## 3)Builed ExpressionSet
library(Biobase)
table(Cluster.df$rownames %in% colnames(COAD.mergedMethyl.sub.M))
sampleID <- Cluster.df$rownames[Cluster.df$rownames %in% colnames(COAD.mergedMethyl.sub.M)]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df[sampleID,])
COAD.methy.set <- ExpressionSet(assayData=COAD.mergedMethyl.sub.M[,sampleID],
                                phenoData=Cluster.df.sub)

















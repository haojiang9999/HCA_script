### 3.GGplot_COAD_miRNA_Batch_Normalized.R
## Using log2(x+1)transformed miRNA value to limma

COAD.pancan.miRNA.mx <- COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset$COAD.pancan.miRNA.mx
COAD.pancan.miRNA.mx <- as.matrix(COAD.pancan.miRNA.mx)
dim(COAD.pancan.miRNA.mx)
hist(COAD.pancan.miRNA.mx)
boxplot(COAD.pancan.miRNA.mx)
summary(COAD.pancan.miRNA.mx[,4])
table(is.na(COAD.pancan.miRNA.mx))
## 1) Log2(x + 1) transformation
COAD.miRNA.log2 <- log2(COAD.pancan.miRNA.mx + 1)
hist(COAD.miRNA.log2)
boxplot(COAD.miRNA.log2)

## 2)Builed ExpressionSet
library(Biobase)
table(Cluster.df$rownames %in% colnames(COAD.miRNA.log2))
sampleID <- Cluster.df$rownames[Cluster.df$rownames %in% colnames(COAD.miRNA.log2)]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df[sampleID,])
COAD.miRNA.set <- ExpressionSet(assayData=COAD.miRNA.log2[,sampleID],
                                phenoData=Cluster.df.sub)
## 3) differential miRNA analysis limma version
## Several Groups experiment design
library(limma)
Cluster.df.sub@data$dynamicColors
f.miRNA <- factor(Cluster.df.sub@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design.miRNA <- model.matrix(~0+f.miRNA)
colnames(design.miRNA) <- c("blue","brown","turquoise","yellow")
# fit the linear model 
fit.miRNA <- lmFit(COAD.miRNA.set, design.miRNA)
contrast.matrix.miRNA <- makeContrasts(blue-(brown+turquoise)/2, 
                                       brown-(blue+turquoise)/2,
                                       turquoise-(brown+blue)/2,
                                       levels=design.miRNA)

# fit the contrasts
fit.miRNA.2 <- contrasts.fit(fit.miRNA, contrast.matrix.miRNA)
fit.miRNA.2 <- eBayes(fit.miRNA.2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.miRNA.2))
# get the table of results for the first contrast blue-(brown+turquoise)
adjPvalueCutoff <- 0.01
logFCcutoff <- log2(2)
number = 50
miRNA.blue <- topTable(fit.miRNA.2, num=number, coef="blue - (brown + turquoise)/2")
miRNA.brown <- topTable(fit.miRNA.2, num=number, coef="brown - (blue + turquoise)/2")
miRNA.turquoise <- topTable(fit.miRNA.2, num=number, coef="turquoise - (brown + blue)/2")
head(miRNA.blue)
## 4) Rescale the data for ploting
COAD.miRNA.set
library(scales)
COAD.miRNA.set.scaled <- t(apply(COAD.miRNA.set, 1, rescale, to=c(-2,2)))
## 5) Heatmap ploting
pheatmap::pheatmap(COAD.miRNA.set.scaled[rownames(miRNA.blue),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "miRNA_Blue VS orthers")
pheatmap::pheatmap(COAD.miRNA.set.scaled[rownames(miRNA.brown),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "miRNA_brown VS orthers")
pheatmap::pheatmap(COAD.miRNA.set.scaled[rownames(miRNA.turquoise),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "miRNA_turquoise VS orthers")

### all samples
pheatmap::pheatmap(COAD.miRNA.set.scaled,annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "miRNA_all")

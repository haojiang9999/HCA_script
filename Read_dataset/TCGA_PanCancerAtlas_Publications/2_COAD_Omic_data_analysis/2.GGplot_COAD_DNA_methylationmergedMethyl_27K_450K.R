### 2.GGplot_COAD_DNA_methylationmergedMethyl_27K_450K.R
# BiocManager::install("lumi")
COAD.mergedMethyl.27.450.mx <- COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset$COAD.mergedMethyl.27.450.mx
dim(COAD.mergedMethyl.27.450.mx)
table(is.na(COAD.mergedMethyl.27.450.mx))
## 1) remove rows contain NAs
COAD.mergedMethyl.sub <- COAD.mergedMethyl.27.450.mx[complete.cases(COAD.mergedMethyl.27.450.mx), ]
dim(COAD.mergedMethyl.sub)
## 2) Convert methylation Beta-value to M-value
library(lumi)
COAD.mergedMethyl.sub.M <- beta2m(COAD.mergedMethyl.sub)
hist(COAD.mergedMethyl.sub.M)
boxplot(COAD.mergedMethyl.sub.M)
class(COAD.mergedMethyl.sub.M)
COAD.mergedMethyl.sub.M[1:5,1:5]
dim(COAD.mergedMethyl.sub.M)
## 3)Builed ExpressionSet
library(Biobase)
table(Cluster.df$rownames %in% colnames(COAD.mergedMethyl.sub.M))
sampleID <- Cluster.df$rownames[Cluster.df$rownames %in% colnames(COAD.mergedMethyl.sub.M)]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df[sampleID,])
COAD.methy.set <- ExpressionSet(assayData=COAD.mergedMethyl.sub.M[,sampleID],
                                phenoData=Cluster.df.sub)

## 4) Probe-wise differential methylation analysis
#https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
## Several Groups experiment design
library(limma)
Cluster.df.sub@data$dynamicColors
f.meth <- factor(Cluster.df.sub@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design.meth <- model.matrix(~0+f.meth)
colnames(design.meth) <- c("blue","brown","turquoise","yellow")
# fit the linear model 
fit.meth <- lmFit(COAD.methy.set, design.meth)
contrast.matrix.methy <- makeContrasts(blue-(brown+turquoise)/2, 
                                      brown-(blue+turquoise)/2,
                                      turquoise-(brown+blue)/2,
                                      levels=design.meth)
# fit the contrasts
fit.meth.2 <- contrasts.fit(fit.meth, contrast.matrix.methy)
fit.meth.2 <- eBayes(fit.meth.2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.meth.2))
# get the table of results for the first contrast blue-(brown+turquoise)
adjPvalueCutoff <- 0.01
logFCcutoff <- log2(2)
number = 50
DMPs.blue <- topTable(fit.meth.2, num=number, coef="blue - (brown + turquoise)/2")
DMPs.brown <- topTable(fit.meth.2, num=number, coef="brown - (blue + turquoise)/2")
DMPs.turquoise <- topTable(fit.meth.2, num=number, coef="turquoise - (brown + blue)/2")
head(DMPs.blue)
## 5) reScale data from -2 to 2 for ploting
library(scales)
COAD.methy.set.scaled <- t(apply(COAD.methy.set, 1, rescale, to=c(-2,2)))

## 6) Heatmap ploting
# Using rescaled data
pheatmap::pheatmap(COAD.methy.set.scaled[rownames(DMPs.blue),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_Blue VS orthers")
pheatmap::pheatmap(COAD.methy.set.scaled[rownames(DMPs.brown),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_brown VS orthers")
pheatmap::pheatmap(COAD.methy.set.scaled[rownames(DMPs.turquoise),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_turquoise VS orthers")



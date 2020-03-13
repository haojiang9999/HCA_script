#### 2.Methylation_data_analysis.R
# It's from TCGA PanCancerAtlas_Publish series paper data
# https://gdc.cancer.gov/about-data/publications/pancanatlas
### 1) Loading DNA methylation mergedMethyl_27K_450K_dataset
COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/DNA_Methylation_Merged_27K_450K_Only/COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset.rds")
COAD.mergedMethyl.27.450.mx <- COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset$COAD.mergedMethyl.27.450.mx
COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset$metadata #
### 2)Read cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
hclust.Res <- Cluster.20200201.V7.Tumor$hclust.Res

## 3) remove rows contain NAs
COAD.mergedMethyl.sub <- COAD.mergedMethyl.27.450.mx[complete.cases(COAD.mergedMethyl.27.450.mx), ]
dim(COAD.mergedMethyl.sub)
## 4) Convert methylation Beta-value to M-value
library(lumi)
COAD.mergedMethyl.sub.M <- beta2m(COAD.mergedMethyl.sub)
hist(COAD.mergedMethyl.sub.M)
boxplot(COAD.mergedMethyl.sub.M)
class(COAD.mergedMethyl.sub.M)
COAD.mergedMethyl.sub.M[1:5,1:5]
dim(COAD.mergedMethyl.sub.M)
## 5)Builed ExpressionSet
library(Biobase)
table(Cluster.df$rownames %in% colnames(COAD.mergedMethyl.sub.M))
sampleID <- Cluster.df$rownames[Cluster.df$rownames %in% colnames(COAD.mergedMethyl.sub.M)]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df[sampleID,])
COAD.methy.set <- ExpressionSet(assayData=COAD.mergedMethyl.sub.M[,sampleID],
                                phenoData=Cluster.df.sub)

## 6) Probe-wise differential methylation analysis
#https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
## Several Groups experiment design
library(limma)
Cluster.df.sub@data$dynamicColors
f.meth <- factor(Cluster.df.sub@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design.meth <- model.matrix(~0+f.meth)
colnames(design.meth) <- c("blue","brown","turquoise","yellow")
## 7) limma DE gene contrast building
fit.meth <- lmFit(COAD.methy.set, design.meth)
contrast.matrix.methy <- makeContrasts(blue-(brown+turquoise)/2, 
                                       brown-(blue+turquoise)/2,
                                       turquoise-(brown+blue)/2,
                                       levels=design.meth)
# fit the contrasts
fit.meth.2 <- contrasts.fit(fit.meth, contrast.matrix.methy)
fit.meth.2 <- eBayes(fit.meth.2)
plotMD(fit.meth.2)
plotMD(fit.meth.2, column = 2)
plotMD(fit.meth.2, column = 3)
## 8) Set up cutoff for top DE genes
adjPvalueCutoff <- 0.01
logFCcutoff <- 0.5
summary(decideTests(fit.meth.2))
number = 100
DMPs.blue <- topTable(fit.meth.2, num=number, coef="blue - (brown + turquoise)/2")
DMPs.brown <- topTable(fit.meth.2, num=number, coef="brown - (blue + turquoise)/2")
DMPs.turquoise <- topTable(fit.meth.2, num=number, coef="turquoise - (brown + blue)/2")
head(DMPs.blue)

## 9) Heatmap ploting
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
pheatmap::pheatmap(COAD.methy.set[rownames(DMPs.blue),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.methy.set[rownames(DMPs.brown),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.methy.set[rownames(DMPs.turquoise),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_turquoise VS orthers",annotation_colors = ann_colors)


### 10) Ploting rescaled methylation data
library(scales)
COAD.methy.scaled <- t(apply(COAD.methy.set, 1, rescale, to=c(-2,2)))
summary(COAD.methy.scaled[5,])
pheatmap::pheatmap(COAD.methy.scaled[rownames(DMPs.blue),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_Blue VS orthers_rescaled",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.methy.scaled[rownames(DMPs.brown),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_brown VS orthers_rescaled",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.methy.scaled[rownames(DMPs.turquoise),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "Methylation_turquoise VS orthers_rescaled",annotation_colors = ann_colors)





















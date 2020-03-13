#### 3.miRNA_data_analysis.R
# It's from TCGA PanCancerAtlas_Publish series paper data
# https://gdc.cancer.gov/about-data/publications/pancanatlas
## Using log2(x+1)transformed miRNA value to limma
### 1) Loading miRNA miRNA_Batch_Normalized_dataset 
COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/miRNA_Batch_Effects_Normalized_miRNA/COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset.rds")
COAD.pancan.miRNA.mx <- COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset$COAD.pancan.miRNA.mx
COAD.pancan.miRNA.mx <- as.matrix(COAD.pancan.miRNA.mx)
COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset$metadata
dim(COAD.pancan.miRNA.mx)
hist(COAD.pancan.miRNA.mx)
boxplot(COAD.pancan.miRNA.mx)
### 2)Read cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
hclust.Res <- Cluster.20200201.V7.Tumor$hclust.Res

### 3) Log2(x + 1) transformation
COAD.miRNA.log2 <- log2(COAD.pancan.miRNA.mx + 1)
hist(COAD.miRNA.log2)
boxplot(COAD.miRNA.log2)

### 4)Builed ExpressionSet
library(Biobase)
table(Cluster.df$rownames %in% colnames(COAD.miRNA.log2))
sampleID <- Cluster.df$rownames[Cluster.df$rownames %in% colnames(COAD.miRNA.log2)]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df[sampleID,])
COAD.miRNA.set <- ExpressionSet(assayData=COAD.miRNA.log2[,sampleID],
                                phenoData=Cluster.df.sub)
## 5) differential miRNA analysis limma version
## Several Groups experiment design
library(limma)
#Cluster.df.sub@data$dynamicColors
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
# Check out ploting
plotMD(fit.miRNA.2)
plotMD(fit.miRNA.2, column = 2)
plotMD(fit.miRNA.2, column = 3)

### 6) Set up cutoff for top DE genes
adjPvalueCutoff <- 0.01
logFCcutoff <- 0.5
summary(decideTests(fit.miRNA.2))
number = 100
miRNA.blue <- topTable(fit.miRNA.2, num=number, coef="blue - (brown + turquoise)/2")
miRNA.brown <- topTable(fit.miRNA.2, num=number, coef="brown - (blue + turquoise)/2")
miRNA.turquoise <- topTable(fit.miRNA.2, num=number, coef="turquoise - (brown + blue)/2")
head(miRNA.blue)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.miRNA.2))

## 7) Heatmap ploting
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
pheatmap::pheatmap(COAD.miRNA.set[rownames(miRNA.blue),],annotation_col = pData(COAD.miRNA.set)[,1:2],
                   main = "miRNA_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.miRNA.set[rownames(miRNA.brown),],annotation_col = pData(COAD.miRNA.set)[,1:2],
                   main = "miRNA_brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.miRNA.set[rownames(miRNA.turquoise),],annotation_col = pData(COAD.miRNA.set)[,1:2],
                   main = "miRNA_turquoise VS orthers",annotation_colors = ann_colors)

### 8) Ploting rescaled miRNAlation data
library(scales)
COAD.miRNA.scaled <- t(apply(COAD.miRNA.set, 1, rescale, to=c(-2,2)))
summary(COAD.miRNA.scaled[5,])
pheatmap::pheatmap(COAD.miRNA.scaled[rownames(miRNA.blue),],annotation_col = pData(COAD.miRNA.set)[,1:2],
                   main = "miRNA_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.miRNA.scaled[rownames(miRNA.brown),],annotation_col = pData(COAD.miRNA.set)[,1:2],
                   main = "miRNA_brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.miRNA.scaled[rownames(miRNA.turquoise),],annotation_col = pData(COAD.miRNA.set)[,1:2],
                   main = "miRNA_turquoise VS orthers",annotation_colors = ann_colors)











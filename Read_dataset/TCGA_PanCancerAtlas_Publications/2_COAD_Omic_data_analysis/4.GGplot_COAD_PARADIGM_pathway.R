### 4.GGplot_COAD_PARADIGM_pathway.R
COAD.PARADIGM.pathway.mx <- COAD_PanCancerAtlas_Publish_PARADIGM_pathway_dataset$COAD.PARADIGM.pathway.mx
dim(COAD.PARADIGM.pathway.mx)
hist(COAD.PARADIGM.pathway.mx)
boxplot(COAD.PARADIGM.pathway.mx)
COAD.PARADIGM.pathway.mx[1:5,1:5]
## Using limma to find significant pathways
## 1)Build an ExpressionSet
library(Biobase)
### convert sample names to pateint barcode
Cluster.df$barcode <- substr(Cluster.df$rownames, start = 1, stop = 12)
table(Cluster.df$barcode %in% colnames(COAD.PARADIGM.pathway.mx))
barcodeID <- Cluster.df$barcode[Cluster.df$barcode %in% colnames(COAD.PARADIGM.pathway.mx)]
rownames(Cluster.df.sub.PARADIGM)
Cluster.df.sub.PARADIGM <- new("AnnotatedDataFrame",
                      data=Cluster.df[Cluster.df$barcode %in% barcodeID,])
rownames(Cluster.df.sub.PARADIGM) <- Cluster.df.sub.PARADIGM@data$barcode
COAD.PARADIGM.set <- ExpressionSet(assayData=COAD.PARADIGM.pathway.mx[,barcodeID],
                                phenoData=Cluster.df.sub.PARADIGM)
## 3) Significant pathway analysis--limma version
## Several Groups experiment design
library(limma)
f.PARADIGM <- factor(Cluster.df.sub.PARADIGM@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design.PARADIGM <- model.matrix(~0+f.PARADIGM)
colnames(design.PARADIGM) <- c("blue","brown","turquoise","yellow")
# fit the linear model 
fit.PARADIGM <- lmFit(COAD.PARADIGM.set, design.PARADIGM)
contrast.matrix.PARADIGM <- makeContrasts(blue-(brown+turquoise), 
                                       brown-(blue+turquoise),
                                       turquoise-(brown+blue),
                                       levels=design.PARADIGM)

# fit the contrasts
fit.PARADIGM.2 <- contrasts.fit(fit.PARADIGM, contrast.matrix.PARADIGM)
fit.PARADIGM.2 <- eBayes(fit.PARADIGM.2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.PARADIGM.2))
# get the table of results for the first contrast blue-(brown+turquoise)
adjPvalueCutoff <- 0.01
logFCcutoff <- log2(2)
number = 50
PARADIGM.blue <- topTable(fit.PARADIGM.2, num=number, coef="blue - (brown + turquoise)")
PARADIGM.brown <- topTable(fit.PARADIGM.2, num=number, coef="brown - (blue + turquoise)")
PARADIGM.turquoise <- topTable(fit.PARADIGM.2, num=number, coef="turquoise - (brown + blue)")
head(PARADIGM.blue)
## 4) Rescale the data for ploting
COAD.PARADIGM.set
library(scales)
COAD.PARADIGM.set.scaled <- t(apply(COAD.PARADIGM.set, 1, rescale, to=c(-2,2)))
## 5) Heatmap ploting
pheatmap::pheatmap(COAD.PARADIGM.set.scaled[rownames(PARADIGM.blue),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "PARADIGM_Blue VS orthers")
pheatmap::pheatmap(COAD.PARADIGM.set.scaled[rownames(PARADIGM.brown),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "PARADIGM_brown VS orthers")
pheatmap::pheatmap(COAD.PARADIGM.set.scaled[rownames(PARADIGM.turquoise),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "PARADIGM_turquoise VS orthers")

pheatmap::pheatmap(COAD.PARADIGM.set[rownames(PARADIGM.blue),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "PARADIGM_Blue VS orthers",scale = "row")
pheatmap::pheatmap(COAD.PARADIGM.set[rownames(PARADIGM.brown),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "PARADIGM_brown VS orthers")
pheatmap::pheatmap(COAD.PARADIGM.set.scaled[rownames(PARADIGM.turquoise),],annotation_col = pData(COAD.methy.set)[,1:2],
                   main = "PARADIGM_turquoise VS orthers")






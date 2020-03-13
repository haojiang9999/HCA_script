#### 2.Significant_Pathway_analysis.R 
### 1.Exame distribution of ssGSEA score
hist(COAD.PanCan33.ssGSEA.xena$TCGA.A6.5657.01)
hist(as.numeric(COAD.PanCan33.ssGSEA.xena[1,]))
COAD.ssGSEA.mx <- as.matrix(COAD.PanCan33.ssGSEA.xena)
plot(density(as.vector(COAD.ssGSEA.mx)), 
     main="COAD.PanCan33.ssGSEA.xena",xlab="ssGSEA score", lwd=2, las=1, xaxt="n", 
     xlim=c(-0.75, 0.75), cex.axis=0.8)
### 2.Using limma to find significant pathway in each groups
library(Biobase)
## 1) Build an ExpressionSet
ssGSEA.score <- as.matrix(COAD.PanCan33.ssGSEA.xena)
head(ssGSEA.score)
summary(ssGSEA.score[,10])
summary(ssGSEA.score[10,])
## How many samples have the ssGSEA score
table(rownames(Cluster.df) %in% colnames(ssGSEA.score))
sampleNames <- rownames(Cluster.df)[rownames(Cluster.df) %in% colnames(ssGSEA.score)]
## Assay data
ssGSEA.score.mx <- ssGSEA.score[,sampleNames]
## Phenotypic data
ssGSEA.pData <- Cluster.df[sampleNames,1:2]
ssGSEA.pData <- new("AnnotatedDataFrame",data=ssGSEA.pData)
table(ssGSEA.pData$dynamicColors)
## Assembling an ExpressionSet
ssGSEASet <- ExpressionSet(assayData=ssGSEA.score.mx,phenoData=ssGSEA.pData)
## Several Groups
f <- factor(ssGSEA.pData@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design <- model.matrix(~0+f)
colnames(design) <- c("blue","brown","turquoise","yellow")
library(limma)
fit <- lmFit(ssGSEASet, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise), 
                                 brown-(blue+turquoise),
                                 turquoise-(brown+blue),
                                   levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH",)
results <- decideTests(fit2)
vennDiagram(results)
results
## Cut off
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(1)

allGeneSets.turquoise <- topTable(fit2, coef="turquoise - (brown + blue)", number=Inf)
DEgeneSets.turquoise <- topTable(fit2, coef="turquoise - (brown + blue)", number=50,
                       p.value=adjPvalueCutoff, lfc=logFCcutoff,adjust="BH")
#res <- decideTests(fit2, p.value=adjPvalueCutoff)
#summary(res)
library(pheatmap)
pheatmap::pheatmap(ssGSEASet[rownames(DEgeneSets.turquoise),],annotation_col = pData(ssGSEASet) )




allGeneSets.blue <- topTable(fit2, coef="blue - (brown + turquoise)", number=Inf)
DEgeneSets.blue <- topTable(fit2, coef="blue - (brown + turquoise)", number=50,
                       p.value=adjPvalueCutoff, lfc=logFCcutoff,adjust="BH")
#res <- decideTests(fit2, p.value=adjPvalueCutoff)
#summary(res)
pheatmap::pheatmap(ssGSEASet[rownames(DEgeneSets.blue),],annotation_col = pData(ssGSEASet) )



allGeneSets.brown <- topTable(fit2, coef="brown - (blue + turquoise)", number=Inf)
DEgeneSets.brown <- topTable(fit2, coef="brown - (blue + turquoise)", number=Inf,
                            p.value=adjPvalueCutoff, lfc=logFCcutoff,adjust="BH")
#res <- decideTests(fit2, p.value=adjPvalueCutoff)
#summary(res)
pheatmap::pheatmap(ssGSEASet[rownames(DEgeneSets.brown),],cluster_cols = hclust.Res, 
                   annotation_col = pData(ssGSEASet) )

##### All ######
# generate some empty ssGSEA samples
y <- exprs(ssGSEASet)
table(hclust.Res$labels %in% colnames(y))
NAsamples <- hclust.Res$labels[!hclust.Res$labels %in% colnames(ssGSEA.score.mx)]
dim(y)
mm <- matrix(0, 1387, 9)
x <- cbind(y,mm)
colnames(x) <- c(colnames(y),NAsamples)
cbind(colnames(y),colnames(x))
NAsamples
pheatmap(x, cluster_cols = hclust.Res)
names(hclust.Res)
hclust.Res$labels

pheatmap(x[rownames(DEgeneSets.blue),], cluster_cols = hclust.Res )
rownames(DEgeneSets.blue)

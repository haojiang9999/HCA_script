##### 4.1_DEgenes_limma_logCPM_methods.R
### 
x <- 2^exprs(COAD.ExpressionSet) - 1
# DGEList object using the edgeR package
library(edgeR)
dge <- DGEList(counts=x)
## filter

keep <- filterByExpr(dge, design.expr)
dge <- dge[keep,,keep.lib.sizes=FALSE]

#normalization to RNA-seq read counts, and the TMM normalization
dge <- calcNormFactors(dge)
dge$samples$norm.factors

# Differential expression: limma-trend
logCPM <- cpm(dge, log=TRUE, prior.count=3)
summary(logCPM)
boxplot(logCPM)
#plotMDS(logCPM)
#
fit.lcpm <- lmFit(logCPM, design.expr)
plotSA(fit.lcpm)
contrast.matrix.expr <- makeContrasts(blue-(brown+turquoise), 
                                      brown-(blue+turquoise),
                                      turquoise-(brown+blue),
                                      levels=design.expr)
fit.lcpm.2 <- contrasts.fit(fit.lcpm, contrast.matrix.expr)
fit.lcpm.2 <- eBayes(fit.lcpm.2)
plotSA(fit.lcpm.2)
plotMD(fit.lcpm.2)

tfit.lcpm <- treat(fit.lcpm.2, lfc=1)
dt.lcpm <- decideTests(tfit.lcpm)
summary(dt.lcpm)
blue.lcpm <- topTreat(tfit.lcpm, coef="blue - (brown + turquoise)", n=Inf)
DEgenes <- rownames(blue.lcpm)
topTreat(tfit.lcpm, coef=1, n=Inf)
basal.vs.ml <- topTreat(dt.lcpm, coef=2, n=Inf)
library(gplots)
heatmap.2(logCPM[DEgenes[1:50],], scale="row",
trace="none", density.info="none",
dendrogram="column")

## 3)Cut off
adjPvalueCutoff <- 0.01
logFCcutoff <- log2(2)
number = 30
DEgenes.lcpm.blue <- topTable(fit.lcpm.2, coef="blue - (brown + turquoise)", number=number,
                         p.value=adjPvalueCutoff,lfc=logFCcutoff,adjust="BH")
pheatmap::pheatmap(logCPM[rownames(DEgenes.lcpm.blue),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_Blue VS orthers",scale="row")



pheatmap::pheatmap(logCPM,annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_Blue VS orthers",scale="row",filename = "logCPM.all.pdf")




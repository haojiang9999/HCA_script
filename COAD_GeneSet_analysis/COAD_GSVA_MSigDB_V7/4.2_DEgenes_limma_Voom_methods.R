### 4.2_DEgenes_limma_Voom_methods.R
x <- 2^exprs(COAD.ExpressionSet) - 1
# DGEList object using the edgeR package
dge <- DGEList(counts=x)
## filter
library(edgeR)
keep <- filterByExpr(dge, design.expr)
dge <- dge[keep,,keep.lib.sizes=FALSE]
## 1) Several Groups experiment design
v <- voom(dge, design.expr,plot=TRUE)
## Normalized Log2 transformed expression data
hist(v$E)
boxplot(v$E)

fit.v <- lmFit(v, design.expr)
contrast.matrix.expr <- makeContrasts(blue-(brown+turquoise)/2, 
                                      brown-(blue+turquoise)/2,
                                      turquoise-(brown+blue)/2,
                                      levels=design.expr)
fit.v.2 <- contrasts.fit(fit.v, contrast.matrix.expr)
fit.v.2 <- eBayes(fit.v.2)
plotSA(fit.v.2)
#tfit <- treat(fit.v.2, lfc=1)
dt <- decideTests(tfit)
summary(dt)
plotMD(tfit, coef="blue - (brown + turquoise)")
## 2) Find top DE genes
adjPvalueCutoff <- 0.01
logFCcutoff <- log2(2)
number = 50
DEgenes.blue.v <- topTable(fit.v.2, coef="blue - (brown + turquoise)", number=number,
                         p.value=adjPvalueCutoff,lfc=logFCcutoff,adjust="BH")
DEgenes.brown.v <- topTable(fit.v.2, coef="brown - (blue + turquoise)", number=number,
                          p.value=adjPvalueCutoff,adjust="BH")
DEgenes.turquoise.v <- topTable(fit.v.2, coef="turquoise - (brown + blue)", number=number,
                              p.value=adjPvalueCutoff,adjust="BH")

## 3) Resacle the expression data
library(scales)
library(scales)
vE.scaled <- t(apply(v$E, 1, rescale, to=c(-2,2)))
summary(vE.scaled[5,])
## 4)Ploting
pheatmap::pheatmap(vE.scaled[rownames(DEgenes.blue.v),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_voom_Blue VS orthers")
pheatmap::pheatmap(vE.scaled[rownames(DEgenes.brown.v),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_voom_brown VS orthers")
pheatmap::pheatmap(vE.scaled[rownames(DEgenes.turquoise.v),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_voom_turquoise VS orthers")
#"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
pheatmap::pheatmap(vE.scaled[rownames(DEgenes.blue.v),],
                   annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_voom_Blue VS orthers", 
                   clustering_distance_rows = "maximum",
                   clustering_distance_cols = "manhattan", 
                   clustering_method = "complete")

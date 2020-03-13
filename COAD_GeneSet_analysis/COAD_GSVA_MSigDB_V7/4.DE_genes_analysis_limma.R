#### 4.DE_genes_analysis_limma.R
### Using limma to find significant DE genes in each groups
library(Biobase)
## 1) Build an ExpressionSet
# COAD.ExpressionSet in 2.GSVA_MSigDB_analysis.R
## 2)Several Groups experiment design
COAD.ExpressionSet@phenoData@data
f.expr <- factor(COAD.ExpressionSet@phenoData@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design.expr <- model.matrix(~0+f.expr)
colnames(design.expr) <- c("blue","brown","turquoise","yellow")

## 4) limma DE analysis
library(limma)
fit.expr <- lmFit(COAD.ExpressionSet, design.expr)
contrast.matrix.expr <- makeContrasts(blue-(brown+turquoise), 
                                 brown-(blue+turquoise),
                                 turquoise-(brown+blue),
                                 levels=design.expr)


fit.expr.2 <- contrasts.fit(fit.expr, contrast.matrix.expr)
fit.expr.2 <- eBayes(fit.expr.2)
## 3)Cut off
adjPvalueCutoff <- 0.01
logFCcutoff <- log2(2)
number = 50
DEgenes.blue <- topTable(fit.expr.2, coef="blue - (brown + turquoise)/2", number=number,
                              p.value=adjPvalueCutoff,lfc=logFCcutoff,adjust="BH")
DEgenes.brown <- topTable(fit.expr.2, coef="brown - (blue + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgenes.turquoise <- topTable(fit.expr.2, coef="turquoise - (brown + blue)/2", number=number,
                                   p.value=adjPvalueCutoff,adjust="BH")
## 5) Heatmap ploting
library(pheatmap)
#pheatmap::pheatmap(COAD.ExpressionSet,annotation_col = pData(COAD.ExpressionSet),
#                   main = "Gene_expression_Blue VS orthers all",show_colnames = F)
pheatmap::pheatmap(COAD.ExpressionSet[rownames(DEgenes.blue),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_Blue VS orthers", scale = "row")
pheatmap::pheatmap(COAD.ExpressionSet[rownames(DEgenes.brown),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_Brown VS orthers", scale = "row")
pheatmap::pheatmap(COAD.ExpressionSet[rownames(DEgenes.turquoise),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_turquoise VS orthers", scale = "row")


cbind(rownames(DEgenes.blue),rownames(DEgenes.brown),rownames(DEgenes.turquoise))


DEgenes.blue <- topTable(fit.expr.2, coef="blue - (brown + turquoise)", number=25,
                         p.value=0.001,lfc=log2(2),adjust="BH")

pheatmap::pheatmap(COAD.ExpressionSet[rownames(DEgenes.blue),],annotation_col = pData(COAD.ExpressionSet),
                   main = "Gene_expression_Blue VS orthers", scale = "row")



y <- voom(COAD.ExpressionSet,design.expr,plot=TRUE)

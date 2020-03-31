#### 3.COAD.Significant_MSigDB_analysis.R
### 1.Exame distribution of ssGSEA score
hist(COAD.c2.all.v7.0)
### 2.Using limma to find significant pathway in each groups
library(Biobase)
## 1) Build an ExpressionSet
COAD.h.all.Set <- ExpressionSet(assayData=COAD.h.all.v7.0,phenoData=Cluster.df.sub)
COAD.c1.all.Set <- ExpressionSet(assayData=COAD.c1.all.v7.0,phenoData=Cluster.df.sub)
COAD.c2.all.Set <- ExpressionSet(assayData=COAD.c2.all.v7.0,phenoData=Cluster.df.sub)
COAD.c3.all.Set <- ExpressionSet(assayData=COAD.c3.all.v7.0,phenoData=Cluster.df.sub)
COAD.c4.all.Set <- ExpressionSet(assayData=COAD.c4.all.v7.0,phenoData=Cluster.df.sub)
COAD.c5.all.Set <- ExpressionSet(assayData=COAD.c5.all.v7.0,phenoData=Cluster.df.sub)
COAD.c6.all.Set <- ExpressionSet(assayData=COAD.c6.all.v7.0,phenoData=Cluster.df.sub)
COAD.c5.all.Set <- ExpressionSet(assayData=COAD.c5.all.v7.0,phenoData=Cluster.df.sub)
## Several Groups experiment design
COAD.c2.all.Set@phenoData@data
f <- factor(Cluster.df.sub@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design <- model.matrix(~0+f)
colnames(design) <- c("blue","brown","turquoise","yellow")
## 2) Significant pathway analysis
###################### COAD.h.all.Set #################
## Cut off
library(limma)
fit.h <- lmFit(COAD.h.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise)/2, 
                                 brown-(blue+turquoise)/2,
                                 turquoise-(brown+blue)/2,
                                 blue-brown,
                                 blue-turquoise,
                                 brown-blue,
                                 brown-turquoise,
                                 turquoise-blue,
                                 turquoise-brown,
                                 levels=design)
fit.h.2 <- contrasts.fit(fit.h, contrast.matrix)
fit.h.2 <- eBayes(fit.h.2)
plotMD(fit.h.2, column = 1)
plotMD(fit.h.2, column = 2)
plotMD(fit.h.2, column = 3)
plotMD(fit.h.2, column = 4)
plotMD(fit.h.2, column = 5)
plotMD(fit.h.2, column = 6)
plotMD(fit.h.2, column = 7)
plotMD(fit.h.2, column = 8)
plotMD(fit.h.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.h.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.h <- topTable(fit.h.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.h <- topTable(fit.h.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.h <- topTable(fit.h.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.h <- topTable(fit.h.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.h <- topTable(fit.h.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.h <- topTable(fit.h.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.h.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "h.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.blue.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.brown.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.turquoise.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.BlvBr.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.BlvTu.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.BrvTu.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_brown - turquoise",annotation_colors = ann_colors)

###################### COAD.c2.all.Set #################
## Cut off
library(limma)
fit.c2 <- lmFit(COAD.c2.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise)/2, 
                                 brown-(blue+turquoise)/2,
                                 turquoise-(brown+blue)/2,
                                 blue-brown,
                                 blue-turquoise,
                                 brown-blue,
                                 brown-turquoise,
                                 turquoise-blue,
                                 turquoise-brown,
                                 levels=design)
fit.c2.2 <- contrasts.fit(fit.c2, contrast.matrix)
fit.c2.2 <- eBayes(fit.c2.2)
plotMD(fit.c2.2, column = 1)
plotMD(fit.c2.2, column = 2)
plotMD(fit.c2.2, column = 3)
plotMD(fit.c2.2, column = 4)
plotMD(fit.c2.2, column = 5)
plotMD(fit.c2.2, column = 6)
plotMD(fit.c2.2, column = 7)
plotMD(fit.c2.2, column = 8)
plotMD(fit.c2.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c2.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c2 <- topTable(fit.c2.2, coef="blue - (brown + turquoise)/2", number=number,
                              p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c2 <- topTable(fit.c2.2, coef="brown - (blue + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c2 <- topTable(fit.c2.2, coef="turquoise - (brown + blue)/2", number=number,
                                   p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c2 <- topTable(fit.c2.2, coef="blue - brown", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c2 <- topTable(fit.c2.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c2 <- topTable(fit.c2.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c2.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c2.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.blue.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.brown.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.turquoise.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.BlvBr.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.BlvTu.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.BrvTu.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_brown - turquoise",annotation_colors = ann_colors)

###################### COAD.c5.all.Set #################
## Cut off
library(limma)
fit.c5 <- lmFit(COAD.c5.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise)/2, 
                                 brown-(blue+turquoise)/2,
                                 turquoise-(brown+blue)/2,
                                 blue-brown,
                                 blue-turquoise,
                                 brown-blue,
                                 brown-turquoise,
                                 turquoise-blue,
                                 turquoise-brown,
                                 levels=design)
fit.c5.2 <- contrasts.fit(fit.c5, contrast.matrix)
fit.c5.2 <- eBayes(fit.c5.2)
plotMD(fit.c5.2, column = 1)
plotMD(fit.c5.2, column = 2)
plotMD(fit.c5.2, column = 3)
plotMD(fit.c5.2, column = 4)
plotMD(fit.c5.2, column = 5)
plotMD(fit.c5.2, column = 6)
plotMD(fit.c5.2, column = 7)
plotMD(fit.c5.2, column = 8)
plotMD(fit.c5.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c5.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c5 <- topTable(fit.c5.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c5 <- topTable(fit.c5.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c5 <- topTable(fit.c5.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c5 <- topTable(fit.c5.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c5 <- topTable(fit.c5.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c5 <- topTable(fit.c5.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c5.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c5.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.blue.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.brown.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.turquoise.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.BlvBr.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.BlvTu.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.BrvTu.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_brown - turquoise",annotation_colors = ann_colors)


###################### COAD.c6.all.Set #################
## Cut off
library(limma)
fit.c6 <- lmFit(COAD.c6.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise)/2, 
                                 brown-(blue+turquoise)/2,
                                 turquoise-(brown+blue)/2,
                                 blue-brown,
                                 blue-turquoise,
                                 brown-blue,
                                 brown-turquoise,
                                 turquoise-blue,
                                 turquoise-brown,
                                 levels=design)
fit.c6.2 <- contrasts.fit(fit.c6, contrast.matrix)
fit.c6.2 <- eBayes(fit.c6.2)
plotMD(fit.c6.2, column = 1)
plotMD(fit.c6.2, column = 2)
plotMD(fit.c6.2, column = 3)
plotMD(fit.c6.2, column = 4)
plotMD(fit.c6.2, column = 5)
plotMD(fit.c6.2, column = 6)
plotMD(fit.c6.2, column = 7)
plotMD(fit.c6.2, column = 8)
plotMD(fit.c6.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c6.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c6 <- topTable(fit.c6.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c6 <- topTable(fit.c6.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c6 <- topTable(fit.c6.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c6 <- topTable(fit.c6.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c6 <- topTable(fit.c6.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c6 <- topTable(fit.c6.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c6.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c6.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.blue.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.brown.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.turquoise.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.BlvBr.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.BlvTu.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.BrvTu.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_brown - turquoise",annotation_colors = ann_colors)



###################### COAD.c1.all.Set #################
## Cut off
library(limma)
fit.c1 <- lmFit(COAD.c1.all.Set, design)
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
fit.c1.2 <- contrasts.fit(fit.c1, contrast.matrix)
fit.c1.2 <- eBayes(fit.c1.2)
plotMD(fit.c1.2, column = 1)
plotMD(fit.c1.2, column = 2)
plotMD(fit.c1.2, column = 3)
plotMD(fit.c1.2, column = 4)
plotMD(fit.c1.2, column = 5)
plotMD(fit.c1.2, column = 6)
plotMD(fit.c1.2, column = 7)
plotMD(fit.c1.2, column = 8)
plotMD(fit.c1.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c1.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c1 <- topTable(fit.c1.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c1 <- topTable(fit.c1.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c1 <- topTable(fit.c1.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c1 <- topTable(fit.c1.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c1 <- topTable(fit.c1.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c1 <- topTable(fit.c1.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c1.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c1.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c1.all.Set[rownames(DEgeneSets.blue.c1),],annotation_col = pData(COAD.c1.all.Set),
                   main = "c1.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c1.all.Set[rownames(DEgeneSets.brown.c1),],annotation_col = pData(COAD.c1.all.Set),
                   main = "c1.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c1.all.Set[rownames(DEgeneSets.turquoise.c1),],annotation_col = pData(COAD.c1.all.Set),
                   main = "c1.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c1.all.Set[rownames(DEgeneSets.BlvBr.c1),],annotation_col = pData(COAD.c1.all.Set),
                   main = "c1.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c1.all.Set[rownames(DEgeneSets.BlvTu.c1),],annotation_col = pData(COAD.c1.all.Set),
                   main = "c1.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c1.all.Set[rownames(DEgeneSets.BrvTu.c1),],annotation_col = pData(COAD.c1.all.Set),
                   main = "c1.all_brown - turquoise",annotation_colors = ann_colors)

###################### COAD.c3.all.Set #################
## Cut off
library(limma)
fit.c3 <- lmFit(COAD.c3.all.Set, design)
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
fit.c3.2 <- contrasts.fit(fit.c3, contrast.matrix)
fit.c3.2 <- eBayes(fit.c3.2)
plotMD(fit.c3.2, column = 1)
plotMD(fit.c3.2, column = 2)
plotMD(fit.c3.2, column = 3)
plotMD(fit.c3.2, column = 4)
plotMD(fit.c3.2, column = 5)
plotMD(fit.c3.2, column = 6)
plotMD(fit.c3.2, column = 7)
plotMD(fit.c3.2, column = 8)
plotMD(fit.c3.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c3.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c3 <- topTable(fit.c3.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c3 <- topTable(fit.c3.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c3 <- topTable(fit.c3.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c3 <- topTable(fit.c3.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c3 <- topTable(fit.c3.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c3 <- topTable(fit.c3.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c3.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c3.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c3.all.Set[rownames(DEgeneSets.blue.c3),],annotation_col = pData(COAD.c3.all.Set),
                   main = "c3.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c3.all.Set[rownames(DEgeneSets.brown.c3),],annotation_col = pData(COAD.c3.all.Set),
                   main = "c3.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c3.all.Set[rownames(DEgeneSets.turquoise.c3),],annotation_col = pData(COAD.c3.all.Set),
                   main = "c3.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c3.all.Set[rownames(DEgeneSets.BlvBr.c3),],annotation_col = pData(COAD.c3.all.Set),
                   main = "c3.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c3.all.Set[rownames(DEgeneSets.BlvTu.c3),],annotation_col = pData(COAD.c3.all.Set),
                   main = "c3.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c3.all.Set[rownames(DEgeneSets.BrvTu.c3),],annotation_col = pData(COAD.c3.all.Set),
                   main = "c3.all_brown - turquoise",annotation_colors = ann_colors)

###################### COAD.c4.all.Set #################
## Cut off
library(limma)
fit.c4 <- lmFit(COAD.c4.all.Set, design)
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
fit.c4.2 <- contrasts.fit(fit.c4, contrast.matrix)
fit.c4.2 <- eBayes(fit.c4.2)
plotMD(fit.c4.2, column = 1)
plotMD(fit.c4.2, column = 2)
plotMD(fit.c4.2, column = 3)
plotMD(fit.c4.2, column = 4)
plotMD(fit.c4.2, column = 5)
plotMD(fit.c4.2, column = 6)
plotMD(fit.c4.2, column = 7)
plotMD(fit.c4.2, column = 8)
plotMD(fit.c4.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c4.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c4 <- topTable(fit.c4.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c4 <- topTable(fit.c4.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c4 <- topTable(fit.c4.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c4 <- topTable(fit.c4.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c4 <- topTable(fit.c4.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c4 <- topTable(fit.c4.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c4.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c4.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c4.all.Set[rownames(DEgeneSets.blue.c4),],annotation_col = pData(COAD.c4.all.Set),
                   main = "c4.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c4.all.Set[rownames(DEgeneSets.brown.c4),],annotation_col = pData(COAD.c4.all.Set),
                   main = "c4.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c4.all.Set[rownames(DEgeneSets.turquoise.c4),],annotation_col = pData(COAD.c4.all.Set),
                   main = "c4.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c4.all.Set[rownames(DEgeneSets.BlvBr.c4),],annotation_col = pData(COAD.c4.all.Set),
                   main = "c4.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c4.all.Set[rownames(DEgeneSets.BlvTu.c4),],annotation_col = pData(COAD.c4.all.Set),
                   main = "c4.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c4.all.Set[rownames(DEgeneSets.BrvTu.c4),],annotation_col = pData(COAD.c4.all.Set),
                   main = "c4.all_brown - turquoise",annotation_colors = ann_colors)

###################### COAD.c7.all.Set #################
## Cut off
library(limma)
fit.c7 <- lmFit(COAD.c7.all.Set, design)
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
fit.c7.2 <- contrasts.fit(fit.c7, contrast.matrix)
fit.c7.2 <- eBayes(fit.c7.2)
plotMD(fit.c7.2, column = 1)
plotMD(fit.c7.2, column = 2)
plotMD(fit.c7.2, column = 3)
plotMD(fit.c7.2, column = 4)
plotMD(fit.c7.2, column = 5)
plotMD(fit.c7.2, column = 6)
plotMD(fit.c7.2, column = 7)
plotMD(fit.c7.2, column = 8)
plotMD(fit.c7.2, column = 9)
## Find Top significant genesets
## Cut off
adjPvalueCutoff <- 0.01
number = 50
## 
dt <- decideTests(fit.c7.2, p.value = adjPvalueCutoff )
summary(dt)
DEgeneSets.blue.c7 <- topTable(fit.c7.2, coef="blue - (brown + turquoise)/2", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c7 <- topTable(fit.c7.2, coef="brown - (blue + turquoise)/2", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c7 <- topTable(fit.c7.2, coef="turquoise - (brown + blue)/2", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvBr.c7 <- topTable(fit.c7.2, coef="blue - brown", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BlvTu.c7 <- topTable(fit.c7.2, coef="blue - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.BrvTu.c7 <- topTable(fit.c7.2, coef="brown - turquoise", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
## Ploting
library(pheatmap)
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
#pheatmap::pheatmap(COAD.c7.all.Set,annotation_col = pData(COAD.h.all.Set),
#                   main = "c7.all_all",show_colnames = F)
pheatmap::pheatmap(COAD.c7.all.Set[rownames(DEgeneSets.blue.c7),],annotation_col = pData(COAD.c7.all.Set),
                   main = "c7.all_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c7.all.Set[rownames(DEgeneSets.brown.c7),],annotation_col = pData(COAD.c7.all.Set),
                   main = "c7.all_Brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c7.all.Set[rownames(DEgeneSets.turquoise.c7),],annotation_col = pData(COAD.c7.all.Set),
                   main = "c7.all_turquoise VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c7.all.Set[rownames(DEgeneSets.BlvBr.c7),],annotation_col = pData(COAD.c7.all.Set),
                   main = "c7.all_Blue - brown",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c7.all.Set[rownames(DEgeneSets.BlvTu.c7),],annotation_col = pData(COAD.c7.all.Set),
                   main = "c7.all_Blue - turquoise",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.c7.all.Set[rownames(DEgeneSets.BrvTu.c7),],annotation_col = pData(COAD.c7.all.Set),
                   main = "c7.all_brown - turquoise",annotation_colors = ann_colors)


#### 3.Significant_MSigDB_analysis.R
### 1.Exame distribution of ssGSEA score
hist(COAD.c2.all.v7.0)
### 2.Using limma to find significant pathway in each groups
library(Biobase)
## 1) Build an ExpressionSet
COAD.h.all.Set <- ExpressionSet(assayData=COAD.h.all.v7.0,phenoData=Cluster.df.sub)
COAD.c2.all.Set <- ExpressionSet(assayData=COAD.c2.all.v7.0,phenoData=Cluster.df.sub)
COAD.c5.all.Set <- ExpressionSet(assayData=COAD.c5.all.v7.0,phenoData=Cluster.df.sub)
COAD.c6.all.Set <- ExpressionSet(assayData=COAD.c6.all.v7.0,phenoData=Cluster.df.sub)
## Several Groups experiment design
COAD.c2.all.Set@phenoData@data
  
f <- factor(Cluster.df.sub@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design <- model.matrix(~0+f)
colnames(design) <- c("blue","brown","turquoise","yellow")
###################### COAD.h.all.Set #################
## Cut off
adjPvalueCutoff <- 0.001
number = 50
library(pheatmap)
library(limma)
fit.h <- lmFit(COAD.h.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise), 
                                 brown-(blue+turquoise),
                                 turquoise-(brown+blue),
                                 levels=design)
fit.h.2 <- contrasts.fit(fit.h, contrast.matrix)
fit.h.2 <- eBayes(fit.h.2)
plotMD(fit.h.2)
DEgeneSets.blue.h <- topTable(fit.h.2, coef="blue - (brown + turquoise)", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.h <- topTable(fit.h.2, coef="brown - (blue + turquoise)", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.h <- topTable(fit.h.2, coef="turquoise - (brown + blue)", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
#library(pheatmap)
pheatmap::pheatmap(COAD.h.all.Set,annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Blue VS orthers all",show_colnames = F)
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.blue.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Blue VS orthers")
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.brown.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_Brown VS orthers")
pheatmap::pheatmap(COAD.h.all.Set[rownames(DEgeneSets.turquoise.h),],annotation_col = pData(COAD.h.all.Set),
                   main = "h.all_turquoise VS orthers")
###################### COAD.c2.all.Set #################
## Cut off
adjPvalueCutoff <- 0.001
number = 50
library(pheatmap)
library(limma)
fit.c2 <- lmFit(COAD.c2.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise), 
                                 brown-(blue+turquoise),
                                 turquoise-(brown+blue),
                                 levels=design)
fit.c2.2 <- contrasts.fit(fit.c2, contrast.matrix)
fit.c2.2 <- eBayes(fit.c2.2)

DEgeneSets.blue.c2 <- topTable(fit.c2.2, coef="blue - (brown + turquoise)", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c2 <- topTable(fit.c2.2, coef="brown - (blue + turquoise)", number=number,
                                 p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c2 <- topTable(fit.c2.2, coef="turquoise - (brown + blue)", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
#library(pheatmap)
pheatmap::pheatmap(COAD.c2.all.Set,annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Blue VS orthers all",show_rownames = F, show_colnames = F)
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.blue.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Blue VS orthers")
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.brown.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_Brown VS orthers")
pheatmap::pheatmap(COAD.c2.all.Set[rownames(DEgeneSets.turquoise.c2),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c2.all_turquoise VS orthers")


###################### COAD.c5.all.Set #################
library(pheatmap)
library(limma)
fit.c5 <- lmFit(COAD.c5.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise), 
                                 brown-(blue+turquoise),
                                 turquoise-(brown+blue),
                                 levels=design)
fit.c5.2 <- contrasts.fit(fit.c5, contrast.matrix)
fit.c5.2 <- eBayes(fit.c5.2)

DEgeneSets.blue.c5 <- topTable(fit.c5.2, coef="blue - (brown + turquoise)", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c5 <- topTable(fit.c5.2, coef="brown - (blue + turquoise)", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c5 <- topTable(fit.c5.2, coef="turquoise - (brown + blue)", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
#library(pheatmap)
## All
pheatmap::pheatmap(COAD.c5.all.Set,annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Blue VS orthers All features",show_rownames = F, show_colnames = F)

pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.blue.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Blue VS orthers")
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.brown.c5),],annotation_col = pData(COAD.c5.all.Set),
                   main = "c5.all_Brown VS orthers")
pheatmap::pheatmap(COAD.c5.all.Set[rownames(DEgeneSets.turquoise.c5),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c5.all_turquoise VS orthers")

###################### COAD.c6.all.Set #################
library(pheatmap)
library(limma)
fit.c6 <- lmFit(COAD.c6.all.Set, design)
contrast.matrix <- makeContrasts(blue-(brown+turquoise), 
                                 brown-(blue+turquoise),
                                 turquoise-(brown+blue),
                                 levels=design)
fit.c6.2 <- contrasts.fit(fit.c6, contrast.matrix)
fit.c6.2 <- eBayes(fit.c6.2)

DEgeneSets.blue.c6 <- topTable(fit.c6.2, coef="blue - (brown + turquoise)", number=number,
                               p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.brown.c6 <- topTable(fit.c6.2, coef="brown - (blue + turquoise)", number=number,
                                p.value=adjPvalueCutoff,adjust="BH")
DEgeneSets.turquoise.c6 <- topTable(fit.c6.2, coef="turquoise - (brown + blue)", number=number,
                                    p.value=adjPvalueCutoff,adjust="BH")
#library(pheatmap)
## All
pheatmap::pheatmap(COAD.c6.all.Set,annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Blue VS orthers All features",show_rownames = F, show_colnames = F)

pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.blue.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Blue VS orthers")
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.brown.c6),],annotation_col = pData(COAD.c6.all.Set),
                   main = "c6.all_Brown VS orthers")
pheatmap::pheatmap(COAD.c6.all.Set[rownames(DEgeneSets.turquoise.c6),],annotation_col = pData(COAD.c2.all.Set),
                   main = "c6.all_turquoise VS orthers")
























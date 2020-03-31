#### 4.RPPA_data_analysis.R
### https://xenabrowser.net/datapages/?dataset=TCGA-RPPA-pancan-clean.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
### 1) Loaging Protein data RPPA
COAD_TCGA_RPPA_clean_xena_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/RPPA_n_7744/COAD_TCGA_RPPA_clean_xena_dataset.rds")
COAD.RPPA.pancan.clean.xena <- COAD_TCGA_RPPA_clean_xena_dataset$COAD.RPPA.pancan.clean.xena
COAD.RPPA.pancan.clean.xena <- as.matrix(COAD.RPPA.pancan.clean.xena)
hist(COAD.RPPA.pancan.clean.xena)
boxplot(COAD.RPPA.pancan.clean.xena)
### 2)Read cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
hclust.Res <- Cluster.20200201.V7.Tumor$hclust.Res

### 3)Builed ExpressionSet
library(Biobase)
table(Cluster.df$rownames %in% colnames(COAD.RPPA.pancan.clean.xena))
sampleID <- Cluster.df$rownames[Cluster.df$rownames %in% colnames(COAD.RPPA.pancan.clean.xena)]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df[sampleID,])
COAD.RPPA.set <- ExpressionSet(assayData=COAD.RPPA.pancan.clean.xena[,sampleID],
                                phenoData=Cluster.df.sub)

## 4) differential RPPA analysis limma version
## Several Groups experiment design
library(limma)
#Cluster.df.sub@data$dynamicColors
f.RPPA <- factor(Cluster.df.sub@data$dynamicColors, levels=c("blue","brown","turquoise","yellow"))
design.RPPA <- model.matrix(~0+f.RPPA)
colnames(design.RPPA) <- c("blue","brown","turquoise","yellow")
# fit the linear model 
fit.RPPA <- lmFit(COAD.RPPA.set, design.RPPA)
contrast.matrix.RPPA <- makeContrasts(blue-(brown+turquoise)/2, 
                                       brown-(blue+turquoise)/2,
                                       turquoise-(brown+blue)/2,
                                       levels=design.RPPA)
# fit the contrasts
fit.RPPA.2 <- contrasts.fit(fit.RPPA, contrast.matrix.RPPA)
fit.RPPA.2 <- eBayes(fit.RPPA.2)
# Check out ploting
plotMD(fit.RPPA.2)
plotMD(fit.RPPA.2, column = 2)
plotMD(fit.RPPA.2, column = 3)
### 5) Set up cutoff for top DE RPPA
adjPvalueCutoff <- 0.01
logFCcutoff <- 0.5
summary(decideTests(fit.RPPA.2))
number = 10
RPPA.blue <- topTable(fit.RPPA.2, num=number, coef="blue - (brown + turquoise)/2")
RPPA.brown <- topTable(fit.RPPA.2, num=number, coef="brown - (blue + turquoise)/2")
RPPA.turquoise <- topTable(fit.RPPA.2, num=number, coef="turquoise - (brown + blue)/2")
head(RPPA.blue)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.RPPA.2))


## 6) Heatmap ploting
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))

pheatmap::pheatmap(COAD.RPPA.set[rownames(RPPA.blue),],annotation_col = pData(COAD.RPPA.set)[,1:2],
                   main = "RPPA_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.RPPA.set[rownames(RPPA.brown),],annotation_col = pData(COAD.RPPA.set)[,1:2],
                   main = "RPPA_brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.RPPA.set[rownames(RPPA.turquoise),],annotation_col = pData(COAD.RPPA.set)[,1:2],
                   main = "RPPA_turquoise VS orthers",annotation_colors = ann_colors)

### 7) Ploting rescaled RPPAlation data
ann_colors = list(dynamicColors = c(blue = "blue",turquoise = "turquoise", 
                                    brown = "brown", yellow = "yellow"))
library(scales)
COAD.RPPA.scaled <- t(apply(COAD.RPPA.set, 1, rescale, to=c(-2,2)))
summary(COAD.RPPA.scaled[5,])
pheatmap::pheatmap(COAD.RPPA.scaled[rownames(RPPA.blue),],annotation_col = pData(COAD.RPPA.set)[,1:2],
                   main = "RPPA_Blue VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.RPPA.scaled[rownames(RPPA.brown),],annotation_col = pData(COAD.RPPA.set)[,1:2],
                   main = "RPPA_brown VS orthers",annotation_colors = ann_colors)
pheatmap::pheatmap(COAD.RPPA.scaled[rownames(RPPA.turquoise),],annotation_col = pData(COAD.RPPA.set)[,1:2],
                   main = "RPPA_turquoise VS orthers",annotation_colors = ann_colors)

### 8) Export tables
RPPA.blue.all <- topTable(fit.RPPA.2, num=Inf, coef="blue - (brown + turquoise)/2")
RPPA.brown.all <- topTable(fit.RPPA.2, num=Inf, coef="brown - (blue + turquoise)/2")
RPPA.turquoise.all <- topTable(fit.RPPA.2, num=Inf, coef="turquoise - (brown + blue)/2")

write.csv(RPPA.blue.all, file = "RPPA.blue.all.csv")
write.csv(RPPA.brown.all, file = "RPPA.brown.all.csv")
write.csv(RPPA.turquoise.all, file = "RPPA.turquoise.all.csv")

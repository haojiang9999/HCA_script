#### 2.Application.R
### 2.Loading data
library(GSEABase)
BiocManager::install("GSVAdata")
library(GSVAdata)
data(c2BroadSets)
c2BroadSets
# We also need to load the following additional libraries:
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
# employ the cache() function from the Biobase package in order to 
#load some pre-computed results and speed up the building time of the vignette:
cacheDir <- system.file("extdata", package="GSVA")
cachePrefix <- "cache4vignette_"
file.remove(paste(cacheDir, list.files(cacheDir, pattern=cachePrefix), sep="/"))

### 2.Functional enrichment
# 1) Filter
data(leukemia)
leukemia_eset
head(pData(leukemia_eset))
table(leukemia_eset$subtype)

filtered_eset <- nsFilter(leukemia_eset, require.entrez=TRUE, remove.dupEntrez=TRUE,
                           var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE,
                           feature.exclude="^AFFX")

leukemia_filtered_eset <- filtered_eset$eset
cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,parallel.sz=10,
                           min.sz=10, max.sz=500, verbose=TRUE),
       dir=cacheDir, prefix=cachePrefix)

leukemia_es
assayData(leukemia_es)
# 2) DE of Gene sets
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)
design <- model.matrix(~ factor(leukemia_es$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgeneSets <- topTable(fit, coef="MLLvsALL", number=Inf,
                          p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)

summary(res)
# When we carry out the corresponding differential expression analysis at gene level:
# 3) DE 0f genes
logFCcutoff <- log2(2)
design <- model.matrix(~ factor(leukemia_eset$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_filtered_eset, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgenes <- topTable(fit, coef="MLLvsALL", number=Inf,
                      p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)
res <- decideTests(fit, p.value=adjPvalueCutoff, lfc=logFCcutoff)
summary(res)

# 4) Plotting
#### DE genessets 
pheatmap::pheatmap(leukemia_es[rownames(DEgeneSets),])









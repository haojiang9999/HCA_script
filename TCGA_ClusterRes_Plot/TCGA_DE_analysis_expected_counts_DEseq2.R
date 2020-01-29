### DE gene cluster resaults
cancerType = "COAD"
COAD.sur.df  # cluster Result
#### 1.Read TCGA dataset ####
#print = TRUE
#geneName = "ENSG00000232956"
filePath <- paste0("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_Pan_TOIL_RSEM_expected_count/",cancerType,"_UCSC_Toil_expected_count_dataset.rds")
TCGA.Toil.expected_count.dataset <- readRDS(filePath)
TCGA.RSEM.gene.expected_count <- TCGA.Toil.expected_count.dataset[[paste0(cancerType,".RSEM.gene.expected_count_round")]]
TCGA.pheno <- TCGA.Toil.expected_count.dataset[[paste0(cancerType,".pheno")]]
TCGA.RSEM.gene.expected_count[1:5,1:5]
#### 2. Select tumor samples
## sample type ###
sampleName <- rownames(TCGA.pheno)
samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
table(samplType)
## "01" was primary tumor, "11' was normal tissue
tumorIndex <- samplType == "01"
#normalIndex <- samplType == "11"
TCGA.pheno.tumor <- TCGA.pheno[tumorIndex,]
tumorSamples <- rownames(TCGA.pheno.tumor)
TCGA.RSEM.gene.expected_count.Tumor <- TCGA.RSEM.gene.expected_count[,tumorSamples]
TCGA.RSEM.gene.expected_count.Tumor[1:5,1:5]
#### 3.DE gene analysis using DEseq2
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(5))
dds <- DESeqDataSetFromMatrix(countData = TCGA.RSEM.gene.expected_count.Tumor,
                              colData = COAD.sur.df[colnames(TCGA.RSEM.gene.expected_count.Tumor),],
                              design= ~ COAD.cluster)
## 1)Pre-filtering
# keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds,parallel = TRUE)
resultsNames(dds) # lists the coefficients
res <- results(dds,name = "COAD.cluster")
summary(res)
res
## 2）Exploring and exporting results
# or to shrink log fold changes association with condition:
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="COAD.cluster", type="apeglm")
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

## 3)p-values and adjusted p-values
# We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]
summary(res)

## 3)Data transformations and visualization
#Extracting transformed values
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
# Effects of transformations on the variance
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
plotDispEsts(dds)
## 4)Data quality assessment by sample clustering and visualization
#Heatmap of the count matrix
library("pheatmap")
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c(1:2)])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,filename = "test.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=FALSE,filename = "test.pdf")
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 500 )
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, cluster_rows=T, show_rownames=F,show_colnames = F,
         cluster_cols=T, annotation_col=df)
#pheatmap(mat, cluster_rows=T, show_rownames=F,
#         cluster_cols=T, annotation_col=df,filename = "test.pdf")


res <- results(vsd,name = "COAD.cluster")

vsdOrdered <- vsd[order(vsd$pvalue),]

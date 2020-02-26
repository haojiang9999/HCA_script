#### 3.Test_on_COAD.R
### 1.Loading Data
TCGA_COAD_RNAseqV2_normalized_log2_dataset <- readRDS("TCGA_COAD_RNAseqV2_normalized_log2_dataset.rds")
COAD.HiSeqV2.log2 <- TCGA_COAD_RNAseqV2_normalized_log2_dataset$COAD.HiSeqV2.log2
COAD.pheno <- TCGA_COAD_RNAseqV2_normalized_log2_dataset$COAD.pheno
### 2.Construct Expression set object
library(Biobase)
all(rownames(COAD.pheno)==colnames(COAD.HiSeqV2.log2))
#colnames(COAD.HiSeqV2.log2)[!colnames(COAD.HiSeqV2.log2) %in% rownames(COAD.pheno)]
#rownames(COAD.pheno)[!rownames(COAD.pheno) %in% colnames(COAD.HiSeqV2.log2)]
## Subset the expression data to COAD.pheno
expr.sub <- as.matrix(COAD.HiSeqV2.log2[,colnames(COAD.HiSeqV2.log2) %in% rownames(COAD.pheno)])
COAD.pheno.sub <- COAD.pheno[colnames(expr.sub),]
COAD.pheno.sub <- new("AnnotatedDataFrame",
                 data=COAD.pheno.sub)
all(rownames(COAD.pheno.sub)==colnames(expr.sub))
COAD.ExpressionSet <- ExpressionSet(assayData=expr.sub,phenoData=COAD.pheno.sub)
### 3.Read GeneSetCollection object gmt download from MSigDB
MSigDB_Geneset_Object_V7_GMT_symbols <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/MSigDB/MSigDB_Geneset_Object_V7_GMT_symbols.rds")


### 4.GSVA
library(GSVA)
COAD.c2.all.v7.0 <- gsva(exprs(COAD.ExpressionSet), MSigDB_Geneset_Object_V7_GMT_symbols$c2.all.v7.0.symbols.gmt,
     parallel.sz=10,min.sz=10, max.sz=500, verbose=TRUE)
COAD.c2.all.v7.0[1:5,1:5]
rownames(COAD.c2.all.v7.0)
table(is.na(COAD.c2.all.v7.0))
dim(COAD.c2.all.v7.0)

### 5.Different gene and pathway analysis
##stringent cut-offs 
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)
# where we will use the latter only for the gene-level differential expression analysis.

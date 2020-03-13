#### 2.COAD_GSVA_MSigDB_analysis.R
## 1)Construct COAD Expression set object
library(Biobase)
all(rownames(Cluster.df)==colnames(COAD.HiSeqV2.log2))
## Subset the expression data to Cluster.df
expr.sub <- as.matrix(COAD.HiSeqV2.log2[,colnames(COAD.HiSeqV2.log2) %in% rownames(Cluster.df)])
Cluster.df.sub <- Cluster.df[colnames(expr.sub),1:2]
Cluster.df.sub <- new("AnnotatedDataFrame",
                      data=Cluster.df.sub)
all(rownames(Cluster.df.sub)==colnames(expr.sub))
COAD.ExpressionSet <- ExpressionSet(assayData=expr.sub,phenoData=Cluster.df.sub)
### 2)GSVA
library(GSVA)
COAD.h.all.v7.0 <- gsva(exprs(COAD.ExpressionSet), MSigDB_gmt_V7_symbols$h.all.v7.0.symbols.gmt,
                        parallel.sz=10,min.sz=10, max.sz=500, verbose=TRUE)
COAD.c2.all.v7.0 <- gsva(exprs(COAD.ExpressionSet), MSigDB_gmt_V7_symbols$c2.all.v7.0.symbols.gmt,
                         parallel.sz=10,min.sz=10, max.sz=500, verbose=TRUE)
dim(COAD.c2.all.v7.0)
COAD.c5.all.v7.0 <- gsva(exprs(COAD.ExpressionSet), MSigDB_gmt_V7_symbols$c5.all.v7.0.symbols.gmt,
                         parallel.sz=10,min.sz=10, max.sz=500, verbose=TRUE)
dim(COAD.c5.all.v7.0)
COAD.c6.all.v7.0 <- gsva(exprs(COAD.ExpressionSet), MSigDB_gmt_V7_symbols$c6.all.v7.0.symbols.gmt,
                         parallel.sz=10,min.sz=10, max.sz=500, verbose=TRUE)
dim(COAD.c6.all.v7.0)

### 3)Save the data
saveRDS(COAD.h.all.v7.0, file = "COAD.GSVA.h.all.v7.0.rds")
saveRDS(COAD.c2.all.v7.0, file = "COAD.GSVA.c2.all.v7.0.rds")
saveRDS(COAD.c5.all.v7.0, file = "COAD.GSVA.c5.all.v7.0.rds")
saveRDS(COAD.c6.all.v7.0, file = "COAD.GSVA.c6.all.v7.0.rds")

# construction a data sets form DEGs analysis
df.TPM
list()
x<- clusterResault[[2]]
y <- x$dynamicColors
names(y) <- rownames(x)
y
list.TPM <- list(tpm = as.matrix(df.TPM), condt = y)
# DEGs function
MAST.TPM.res <- run_MASTtpm(list.TPM)
x <- MAST.TPM.res$res
class(x)
x[,"hurdle",]
head(MAST.TPM.res$df)
head(MAST.TPM.res$res)
class(MAST.TPM.res$res)
MAST.TPM.res[,]
x <- MAST.TPM.res$res
x$pri 
MAST.TPM.res$res[, "hurdle"]
# find DEGs by p-values
head(x)
x <- MAST.TPM.res$df
y <- x[order(x$pval),, drop = FALSE]
head(y,100)
class(x)




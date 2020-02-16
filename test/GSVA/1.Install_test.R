#### 1.Install_test.R
BiocManager::install("GSVA")

#### Test
library(GSVA)
p <- 20000 ## number of genes
n <- 30 ## number of samples
nGS <- 100 ## number of gene sets
min.sz <- 10 ## minimum gene set size
max.sz <- 100 ## maximum gene set size
X <- matrix(rnorm(p*n), nrow=p, dimnames=list(1:p, 1:n))
dim(X)

gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE)) ## sample gene set sizes
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets

es.max <- gsva(X, gs, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif <- gsva(X, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
par(mfrow=c(1,2), mar=c(4, 4, 4, 1))
plot(density(as.vector(es.max)), main="Maximum deviation from zero",
        xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",
        xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
       axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
















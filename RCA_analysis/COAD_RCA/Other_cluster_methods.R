## plot using other methods
## plot using correlation directly
set.seed(123)
library(Rtsne)
t(test_for_clust_power)
TSNE <- Rtsne(t(test_for_clust_power), dims = 2, initial_dims = 50, perplexity = 50)
tSNEdata.used <- as.matrix(TSNE$Y)
plot(tSNEdata.used[,1],tSNEdata.used[,2]) ## OK
## filter
test_for_clust_power


######## observation
PCA <- prcomp(t(test_for_clust_power))
PCAdata <- as.matrix(PCA$x[,1:44])
plot(1:20,PCA$sdev[1:20])


## plot
library(ggplot2)
theme_set(theme_bw())
p0419 <- data.frame(tSNE1 = tSNEdata.used[,1], tSNE2 = tSNEdata.used[,2], 
                    Cluster = factor(cluster_result), 
                    PC1 = PCAdata[,1], PC2 = PCAdata[,2], PC3 = PCAdata[,3], PC4 = PCAdata[,4])
plot0419 <- ggplot(p0419,aes(tSNE1,tSNE2,colour=Cluster)) + geom_point()
print(plot0419)

#### Merge data analysis
source("/data8t_4/JH/MyJobs/1_R_script/RCA/RCA_seperate_fuction.R")
test1 <- RCA.cluster(Cor.merged)
test2 <- RCA.PCA(Cor.merged, test1[[2]]$dynamicColors)
test3 <- RCA.heatmap(Cor.merged,test1[[1]], test1[[2]]$dynamicColors)

#### Multiple plot of the data
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/DimeReduPlot.R")
DimeReduPlot(mx = scale(Cor.merged), color = Pheno.merged$Cell_info, tiltle = "No transform")


#### PCA ggplot2
Cor.merged.transform = abs(Cor.merged)^power * sign(Cor.merged)


pca <- prcomp(t(scale(Cor.merged))); 

pca$rotation # feartures contribution to each PCs
pca$x        # samples position on that direction PCs
pca$sdev     # sd of samples on the PCs direction
pca_out <- as.data.frame(pca$x)
#row.names(pca_out) <- gsub("_scTrioSeq2Rna_scTrioSeq2Rna_", "_scTrioSeq2Rna_", row.names(pca_out))
#pca_out$cellType <- sapply(strsplit(as.character(row.names(pca_out)), "_"), "[[", 4 )
#pca_out$pateint <- sapply(strsplit(as.character(row.names(pca_out)), "_"), "[[", 3 )

library(ggplot2)
p<-ggplot(pca_out,aes(x=PC1,y=PC2,color = Pheno.merged$Cell_info ))
p<-p+geom_point()
p
pc_projection = as.data.frame(pr$x);
cell_projection = pc_projection[,1:2];

### t-SNE plot
library(Rtsne)
#t(test_for_clust_power)
TSNE <- Rtsne(t(Cor.merged), dims = 2, initial_dims = 50, perplexity = 50)
tSNEdata.used <- as.matrix(TSNE$Y)
library(ggplot2)
tsne_plot <- data.frame(x = tSNEdata.used[,1], y =tSNEdata.used[,2], col = Pheno.merged$Cell_info)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + labs(title = "t-SNE plot")


#### UMAP plot
library(umap)
UMAP <- umap(t(Cor.merged))
UMAP.out <- UMAP$layout
umap_plot <- data.frame(x =UMAP.out[,1] , y =UMAP.out[,2], col = Pheno.merged$Cell_info)
ggplot(umap_plot) + geom_point(aes(x=x, y=y, color=col))+labs(title = "UMAP plot")

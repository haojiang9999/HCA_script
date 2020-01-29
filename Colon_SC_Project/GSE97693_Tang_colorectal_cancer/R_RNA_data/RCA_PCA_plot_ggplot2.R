#### PCA ggplot2
pca <- prcomp(t(scale(test_for_clust_power))); 

pca$rotation # feartures contribution to each PCs
pca$x        # samples position on that direction PCs
pca$sdev     # sd of samples on the PCs direction
pca_out <- as.data.frame(pca$x)
row.names(pca_out) <- gsub("_scTrioSeq2Rna_scTrioSeq2Rna_", "_scTrioSeq2Rna_", row.names(pca_out))
pca_out$cellType <- sapply(strsplit(as.character(row.names(pca_out)), "_"), "[[", 4 )
pca_out$cellType2 <-substr(pca_out$cellType, start = 1, stop = 2)
pca_out$pateint <- sapply(strsplit(as.character(row.names(pca_out)), "_"), "[[", 3 )
saveRDS(pca_out, file = "pca_TPM_out_table.rds")
library(ggplot2)
p<-ggplot(pca_out,aes(x=PC1,y=PC2,color=pateint ))
p<-p+geom_point()
p
library(ggplot2)
p<-ggplot(pca_out,aes(x=PC2,y=PC3,color = clusterResault[[2]]$dynamicColors ))
p<-p+geom_point()
p

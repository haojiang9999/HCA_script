#### Using RCA pipline for COAD
### calculate the correlation matrix
#test <- corMatrix(log10(df.FPKM+1),refPanel.filtered_of_fetal, method = "pearson")
#test <- corMatrix(GSE97693_Tang_TPM_cells,referPanel_of_mix, method = "pearson")
#saveRDS(test, file = "corMatrix_of_df.TPM.rds")
# convert the cor value
COAD_AnnoRes_Combine <- Combine_df
rownames(COAD_AnnoRes_Combine)
test_for_clust_power = abs(COAD_AnnoRes_Combine)^4 * sign(COAD_AnnoRes_Combine)
#saveRDS(test_for_clust_power, file = "test_for_clust_power.rds")
#getwd()
### Generate cell clusters
clusterResault <- RCA.cluster(test_for_clust_power, deepSplit_wgcna=0, min_group_Size_wgcna=1)
table(clusterResault[[2]]$dynamicColors)
#clusterResault[[1]]
#clusterResault[[2]]
#saveRDS(clusterResault[[2]], file = "clusterResault.TPM.rds")
## add sample imfo
### Function to plot PCA figure
RCA.PCA(test_for_clust_power, color = clusterResault[[2]]$dynamicColors)

### Function to plot heatmap figure
RCA.heatmap(test_for_clust_power, cellTree = clusterResault[[1]], 
            color = clusterResault[[2]]$dynamicColors)

heatmap.2(as.matrix(test_for_clust_power),
          col=color_scheme,
          Colv=as.dendrogram(clusterResault[[1]]),
          ColSideColors=as.vector(clusterResault[[2]]$dynamicColors),
          scale= "column", 
          margins=c(5,28),
          trace="none",
          key = TRUE,
          keysize = 1.5,
          cexCol = 1,cexRow =0.95,
          labCol = "",
          density.info=c("none"))

dev.off()

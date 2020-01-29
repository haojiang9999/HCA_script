#### Using RCA pipline
### calculate the correlation matrix
test <- corMatrix(log10(df.FPKM+1),refPanel.filtered_of_fetal, method = "pearson")
test <- corMatrix(df.TPM,refPanel.filtered_of_fetal, method = "pearson")
saveRDS(test, file = "corMatrix_of_df.TPM.rds")
# convert the cor value
test_for_clust_power = abs(test)^4 * sign(test)
#saveRDS(test_for_clust_power, file = "test_for_clust_power.rds")
#getwd()
### Generate cell clusters
clusterResault <- RCA.cluster(test_for_clust_power, deepSplit_wgcna=0.5, min_group_Size_wgcna=1)
table(clusterResault[[2]]$dynamicColors)
#clusterResault[[1]]
#clusterResault[[2]]
saveRDS(clusterResault[[2]], file = "clusterResault.TPM.rds")
## add sample imfo
### Function to plot PCA figure
RCA.PCA(test_for_clust_power, color = clusterResault[[2]]$dynamicColors)

### Function to plot heatmap figure
RCA.heatmap(test_for_clust_power, cellTree = clusterResault[[1]], 
            color = clusterResault[[2]]$dynamicColors)
dev.off()

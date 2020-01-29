#### Expression as input
# using TPM as input
TPM.fil <- df.TPM[rowSums(df.TPM) > 3,]
# find most varince genes 
CV <- function(x){
  (sd(x)/mean(x))*100
}
genesSelected <- apply(TPM.fil, 1, CV)
genesSelected <- tail(sort(genesSelected), 8000)
TPM.fil <- TPM.fil[names(genesSelected), ]

### Generate cell clusters
clusterResault <- RCA.cluster(TPM.fil)
clusterResault[[1]]

### Function to plot PCA figure
RCA.PCA(TPM.fil, color = clusterResault[[2]]$dynamicColors)
### Function to plot heatmap figure
#RCA.heatmap(TPM.fil, cellTree = clusterResault[[1]], color = clusterResault[[2]]$dynamicColors)
dev.off()
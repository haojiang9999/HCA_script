#### AnnoRes.688 cluster analysis
Combine_df[,1]
# convert the cor value
test_for_clust_power = abs(Combine_df)^4 * sign(Combine_df)
test_for_clust_power[,1]


### PCA plot
#install.packages("ggfortify")
library("ggfortify")

Pheno.688.df <- data.frame(cellType = Pheno.688$cellType,
                           cellType2 = Pheno.688$cellType2,
                           pateint = Pheno.688$pateint,
                           sampleNames = Pheno.688$sampleNames)
colour = "cellType2"

autoplot(prcomp(t(test_for_clust_power)), data = Pheno.688.df, colour = colour )
### T-SNE plot
library(Rtsne)
library(ggplot2)
set.seed(123)
tsne_out <- Rtsne(t(Combine_df), dims = 2, initial_dims = 50, perplexity = 50)
col = Pheno.688$cellType2
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], 
                        col = col # change phenotypes
)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))


#### cluster 
source("/data8t_4/JH/MyJobs/RCA_analysis/Rs_colon_cancer_RCA/RCA_seperate_fuction.R")
clusterResault <- RCA.cluster(test_for_clust_power, deepSplit_wgcna=0, min_group_Size_wgcna=1)
table(clusterResault[[2]]$dynamicColors)
RCA.PCA(test_for_clust_power, color = clusterResault[[2]]$dynamicColors)
RCA.heatmap(test_for_clust_power, cellTree = clusterResault[[1]], 
            color = clusterResault[[2]]$dynamicColors)

#### Cluster using SC3
### Step1.Load script
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/JH_SC3_cluster.R")
### Step2.Convert NA in matrix to 0
Cor.Norm <- Min_max_norm.1500
Cor.Norm[is.na(Cor.Norm)] <-0
colnames(Cor.Norm)[colSums(is.na(Cor.Norm)) > 0]
### Step3.Cluster cells (In my script the cells with 0 variance was removed)
SC3_cluster_1500 <- JH_SC3_cluster(Cor.Norm,Pheno.merged,ks=2:4)
## Select a cluster results
hc_1500 <- SC3_cluster_1500$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Min_max_norm.1500,show_colnames = F,
           annotation_col = Pheno.merged[,c(1,2)], cluster_cols = test2)
### Step4.Group the clusters
## cut at height == 100
test <- as.dendrogram(hc_1500)
test2 <- as.hclust(test)
cutree(hc_1500, 2)
d <- cut(as.dendrogram(hc_1500), h=15)
plot(d[["upper"]])
plot(d[["lower"]][[1]])
plot(d[["lower"]][[2]])

d1 <- d$lower[[2]]
### Extract label in the dendrogram
d1_data <- dendro_data(d1)
head(d1_data$segments)
head(d1_data$labels)

# load package
detach("dendextend")
library(dendextend)
as.dendrogram(hc_1500) %>%
  plot(horiz = F)
ClusterGroup <- cutree(tree = as.dendrogram(hc_1500), k = 2)
table(ClusterGroup)
group_1 <- names(ClusterGroup[ClusterGroup == 1])
group_2 <- names(ClusterGroup[ClusterGroup == 2])
Cor_res_sub_1 <- Cor.Norm[,group_1]
Pheno.merged_1 <- Pheno.merged[group_1,]
heatmap.JH(Cor_res_sub_1,show_colnames = F,
           annotation_col = Pheno.merged_1, cluster_cols = d1)
library(heatmap)

Cor_res_sub_2 <- Cor.Norm[,group_2]
Pheno.merged_2 <- Pheno.merged[group_2,]
heatmap.JH(Cor_res_sub_2,show_colnames = F,
           annotation_col = Pheno.merged_2)

Min_max_norm.1500[,group_x]
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Min_max_norm.1500,show_colnames = F,
           annotation_col = Pheno.merged[,c(1,2)], cluster_cols = hc_1500)

hc_1500


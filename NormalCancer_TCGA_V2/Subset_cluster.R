### subset cluster
Cor.tumor.Adult <- Cor.tumor[1:8,]
Cluster.Adult <- pheatmap::pheatmap(Cor.tumor.Adult,show_colnames = T,cluster_rows = F,scale = "none")
source("/data8t_4/JH/MyJobs/1_R_script/TCGA_plot/TCGAClusterSurv.R")
TCGAClusterSurv(Input.tb = Cor.tumor.Adult, 
                hclust.res = Cluster.Adult$tree_col, Col.anno = Col.anno, k = 2)
# SC3 cluster

Cluster.Adult.SC3 <- JH_SC3_cluster(Cor.tumor.Adult,Pheno.merged.tumor,ks=2:4)
hcTumor <- SC3.Tumor$`3`$hc
TCGAClusterSurv(Input.tb = Cor.tumor.Adult, 
                hclust.res = Cluster.Adult.SC3$`3`$hc, Col.anno = Col.anno, k = 2)

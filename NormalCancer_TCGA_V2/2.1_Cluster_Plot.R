#### 2.1_Cluster_Plot
### Step1.Figure out the most appropriate dissimilarity measure 
hclust.Res <- hclust(dist(t(Trans.Rang1.cv2500_colon_InNormal), method="euclidean"), method="complete")
## Not good for dist(): "maximum","canberra","binary"
## Good for dist():"euclidean"="minkowski">"manhattan"
options(repr.plot.width=20, repr.plot.height=10)
Cluster.Res <- pheatmap::pheatmap(Trans.Rang1.cv2500_colon_InNormal,annotation_col = COAD.pheno.tumor,
                                  cluster_cols = hclust.Res,
                                  show_colnames = F,cluster_rows = F,scale = "none")
source("/data8t_4/JH/MyJobs/1_R_script/TCGA_plot/TCGAClusterSurv.R")
x <- TCGAClusterSurv(Input.tb = Trans.Rang1.cv2500_colon_InNormal, hclust.res = hclust.Res, Col.anno = COAD.pheno[TumorID,], k = 4)
### Cut tree resault
cutree.res <- cutree(hclust.Res, k = 4)
table(cutree.res)
library(WGCNA)
dynamicColors = labels2colors(cutree.res) # convert label to color
Cluster.20200201.V7.Tumor <- list(cutree.res = cutree.res,
                                  dynamicColors = dynamicColors,
                                  hclust.Res = hclust.Res)
saveRDS(Cluster.20200201.V7.Tumor, file = "Cluster.20200201.V7.Tumor.rds")

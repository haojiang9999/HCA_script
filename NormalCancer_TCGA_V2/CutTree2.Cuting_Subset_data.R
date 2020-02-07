# CutTree2.Cuting_test_data.R
### Subset the Trans.Rang1.cv2500 use only the colon info
rownames(Trans.Rang1.cv2500)
Trans.Rang1.cv2500_colon <- Trans.Rang1.cv2500[c(1:39),]
rownames(Trans.Rang1.cv2500_colon)
### Step1.Figure out the most appropriate dissimilarity measure 
hclust.Res <- hclust(dist(t(Trans.Rang1.cv2500_colon), method="euclidean"), method="complete")
## Not good for dist(): "maximum","canberra","binary"
## Good for dist():"euclidean"="minkowski">"manhattan"

Cluster.Res <- pheatmap::pheatmap(Trans.Rang1.cv2500_colon,annotation_col = COAD.pheno.tumor,
                                  cluster_cols = hclust.Res,
                                  show_colnames = F,cluster_rows = F,scale = "none")
#### Step2.Plot dendogram
plot(hclust.Res)

#### Step3.Use the hybrid tree cut method and tune shape parameters
library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(reshape2)
dynamicCut <- cutreeDynamic(hclust.Res, minClusterSize=20, method="tree",
                            distM=as.matrix(dist(t(Trans.Rang1.cv2500)), 
                                            deepSplit=1, maxCoreScatter=NULL, minGap=NULL, 
                                            maxAbsCoreScatter=NULL, minAbsGap=NULL))
#dynamicCut
# Plot dendrogram
plotDendroAndColors(dendro=hclust.Res,colors=dynamicCut)

#### Step4.Use cut tree 
test1<- cutree(hclust.Res, k=3)
test1
source("/data8t_4/JH/MyJobs/1_R_script/TCGA_plot/TCGAClusterSurv.R")
TCGAClusterSurv(Input.tb = Trans.Rang1.cv2500, hclust.res = hclust.Res, Col.anno = COAD.pheno[TumorID,], k = 4)


TCGAClusterSurv.2(Input.tb = Trans.Rang1.cv2500,hclust.res= hclust.Res, group.res = dynamicCut, Col.anno = COAD.pheno[TumorID,])

TCGAClusterSurv.2

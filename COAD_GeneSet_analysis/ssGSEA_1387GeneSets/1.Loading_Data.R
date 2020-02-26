#### 1.Loading Data 
### 1)Loading cluster results
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
### 2)Loading ssGSEA download from Xena
COAD_ssGSEA_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_panCancer_pathway_activity/COAD_TCGA_PanCan33_ssGSEA_dataset.rds")
COAD.PanCan33.ssGSEA.xena <- COAD_ssGSEA_dataset$COAD.PanCan33.ssGSEA.xena
class(COAD.PanCan33.ssGSEA.xena)








#### 1.Loading_data.R
### 1)Read cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
hclust.Res <- Cluster.20200201.V7.Tumor$hclust.Res

### 2) DNA methylation mergedMethyl_27K_450K_dataset
COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/DNA_Methylation_Merged_27K_450K_Only/COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset.rds")
COAD_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset$metadata

### 3) miRNA miRNA_Batch_Normalized_dataset 
COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/miRNA_Batch_Effects_Normalized_miRNA/COAD_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset.rds")

### 4) PARADIGM_pathway_dataset
COAD_PanCancerAtlas_Publish_PARADIGM_pathway_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/PARADIGM_Pathway_Inference_Matrix/COAD_PanCancerAtlas_Publish_PARADIGM_pathway_dataset.rds")

### 5) RNA_Final_dataset
COAD_PanCancerAtlas_Publish_RNA_Final_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/RNA_Final/COAD_PanCancerAtlas_Publish_RNA_Final_dataset.rds")

### 6) RPPA_Final_dataset
COAD_PanCancerAtlas_Publish_RPPA_Final_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/RPPA_Final/COAD_PanCancerAtlas_Publish_RPPA_Final_dataset.rds")



























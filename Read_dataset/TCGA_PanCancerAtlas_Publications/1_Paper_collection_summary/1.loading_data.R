#### 1.loading_data.R 
### 1)Read cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
hclust.Res <- Cluster.20200201.V7.Tumor$hclust.Res
### 2)GI_Adenocarcinomas_Characteristics_And_Epigenetic_Silencing
## Comparative_Molecular_Analysis_of_Gastrointestinal_Adenocarcinomas
file1 <- "/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/Comparative_Molecular_Analysis_of_Gastrointestinal_Adenocarcinomas/COAD_GI_Adenocarcinomas_Characteristics_And_Epigenetic_Silencing_dataset.rds"
COAD_GI_Adenocarcinomas_dataset <- readRDS(file1)
COAD.GI.Adenocarcinomas.Characteristics <- COAD_GI_Adenocarcinomas_dataset$COAD.GI.Adenocarcinomas.Characteristics
COAD.Epigenetic.Silencing.Calls <- COAD_GI_Adenocarcinomas_dataset$COAD.Epigenetic.Silencing.Calls
### 3)COAD_Aneuploidy_score_dataset
## Genomic_and_Functional_Approaches_to_Understanding_Cancer_Aneuploidy
file2 <- "/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/Genomic_and_Functional_Approaches_to_Understanding_Cancer_Aneuploidy/COAD_Aneuploidy_score_dataset.rds"
COAD_Aneuploidy_score_dataset <- readRDS(file2)
COAD_Aneuploidy_score_dataset$COAD.Aneuploidy.score
### 4)COAD_Oncogenic_Signaling_Alteration_dataset
##Oncogenic_Signaling_Pathways_in_The_Cancer
file3 <- "/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/Oncogenic_Signaling_Pathways_in_The_Cancer/COAD_Oncogenic_Signaling_Alteration_dataset.rds"
COAD_Oncogenic_Signaling_Alteration_dataset <- readRDS(file3)
COAD_Oncogenic_Signaling_Alteration_dataset$COAD.Alteration.pathway.level
### Is these NA mean the was deleted?
COAD_Oncogenic_Signaling_Alteration_dataset$COAD.Alteration.gene.level$CDKN2A

### 5)COAD_The_Immune_Landscape_of_Cancer_Immune_Characteristics
##The_Immune_Landscape_of_Cancer
file5 <- "/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/The_Immune_Landscape_of_Cancer/COAD_The_Immune_Landscape_of_Cancer_Immune_Characteristics.rds"
COAD_Immune_Characteristics <- readRDS(file5)























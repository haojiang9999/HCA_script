#### 6.Xena_panCancer_signature_analysis.R
### 1)Loading cluster resaults
Cluster.20200201.V7.Tumor <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_TCGA_V2/Cluster.20200201.V7.Tumor.rds")
cutree.res <- Cluster.20200201.V7.Tumor$cutree.res
dynamicColors <- Cluster.20200201.V7.Tumor$dynamicColors
Cluster.df <- cbind(cutree.res,dynamicColors) 
Cluster.df <- as.data.frame(Cluster.df)
Cluster.df$rownames <- rownames(Cluster.df)
############ Stemness score analysis ############
### 2)Loading StemnessScore dataset
## Paper : Machine Learning Identifies Stemness Features
Associated with Oncogenic Dedifferentiation
## DNAss
COAD_StemnessScores_DNAmeth_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/signatures/Stemness_score_DNA_methylation_based/COAD_StemnessScores_DNAmeth_dataset.rds")
COAD.StemnessScores.DNAmeth<- COAD_StemnessScores_DNAmeth_dataset$COAD.pancancer.StemnessScores.DNAmeth.xena
COAD.StemnessScores.DNAmeth <- as.data.frame(t(COAD.StemnessScores.DNAmeth))
COAD.StemnessScores.DNAmeth$rownames <- rownames(COAD.StemnessScores.DNAmeth)
## RNAss
COAD_StemnessScores_RNAexp_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/signatures/Stemness_score_RNA_based/COAD_StemnessScores_RNAexp_dataset.rds")
COAD.StemnessScores.RNAexp<- COAD_StemnessScores_RNAexp_dataset$COAD.pancancer.StemnessScores.RNAexp.xena
COAD.StemnessScores.RNAexp <- as.data.frame(t(COAD.StemnessScores.RNAexp))
COAD.StemnessScores.RNAexp$rownames <- rownames(COAD.StemnessScores.RNAexp)

########### Immune subtype analysis ########################
## Paper:The Immune Landscape of Cancer
COAD_Subtype_Immune_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_panCancer_phenotype_Immune_subtype/COAD_Subtype_Immune_dataset.rds")
COAD.Subtype.Immune <- COAD_Subtype_Immune_dataset$COAD.Subtype.Immune.xena

### 2)Constucting Ploting dataframe

Plot.df.TCGA.signature <- dplyr::left_join(Cluster.df, COAD.StemnessScores.DNAmeth, by = "rownames")
Plot.df.TCGA.signature <- dplyr::left_join(Plot.df.TCGA.signature, COAD.StemnessScores.RNAexp, by = "rownames")
Plot.df.TCGA.signature <- dplyr::left_join(Plot.df.TCGA.signature, COAD.Subtype.Immune, by = "rownames")

### 3)Plotting
## Stemness DNA RNA
Plot.df.TCGA.signature$RNAss
library(reshape2)
dat.m.stemness <- melt(Plot.df.TCGA.signature[,c(2,4:9)],id.vars='dynamicColors')
ggplot(dat.m.stemness) + geom_boxplot(aes(x=variable, y=value, color=dynamicColors)) +
  scale_color_manual(values= c("blue","brown","turquoise","yellow"))
##Immune subtype
library(reshape2)
Plot.df.TCGA.signature$Subtype_Immune_Model_Based
ggplot(data=subset(Plot.df.TCGA.signature, !is.na(Subtype_Immune_Model_Based)), aes(x = dynamicColors, fill = Subtype_Immune_Model_Based)) + 
  geom_bar(position = "fill") + theme_minimal()+
  scale_fill_manual(values= alpha(c("#999999", "#E69F00", "#56B4E9","#76EE00","#9A32CD","#CD6600"),0.6)) # transparency










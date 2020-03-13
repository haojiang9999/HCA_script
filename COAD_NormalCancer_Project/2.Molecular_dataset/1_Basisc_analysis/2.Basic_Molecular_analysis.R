#### 2.Basic_Molecular_analysis.R
### 1)Total Mutation Plot for each group
#Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines
# Loadingdata
COAD_mc3.v0.2.8.PUBLIC.maf_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/Mutations_mc3_v0.2.8_PUBLIC/COAD_mc3.v0.2.8.PUBLIC.maf_dataset.rds")
COAD_mc3.v0.2.8.PUBLIC.maf_dataset$mc3.v0.2.8.PUBLIC.metadata
COAD.mc3.PUBLIC_sampleSummary <- COAD_mc3.v0.2.8.PUBLIC.maf_dataset$COAD.mc3.PUBLIC_sampleSummary
COAD.mc3.PUBLIC_sampleSummary$rownames <- rownames(COAD.mc3.PUBLIC_sampleSummary)
# Merge table
MergeTable.totalMut <- dplyr::left_join(Cluster.df, COAD.mc3.PUBLIC_sampleSummary, by = "rownames")
# Plot 
library(ggplot2)
# Remove NA value
data=subset(MergeTable.totalMut, !is.na(MergeTable.totalMut$total))
xlabs <- paste(levels(data$dynamicColors),"\n(N=",table(data$dynamicColors),")",sep="")
ggplot(data,aes(x = dynamicColors, y = total,fill=dynamicColors))+
  geom_boxplot()+scale_x_discrete(labels=xlabs) + theme_minimal()+
  scale_fill_manual(values=c("blue","brown","turquoise","yellow")) +
  labs(title =" MC3 totalMut")



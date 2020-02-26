#### 1.Generate_The_Immune_Landscape_of_Cancer_Immune_Characteristics.R
# Paper:The Immune Landscape of Cancer
Immune_Characteristics <- read.csv("Table S1. PanImmune Feature Matrix of Immune Characteristics. Related to Figure 1.csv")
table(Immune_Characteristics$TCGA.Study) ## each cancer types
TCGA.cancer.types <- names(table(Immune_Characteristics$TCGA.Study) )
for(i in TCGA.cancer.types){
  # i= "ACC"
  ### Step1 separate data by cancer types ###
  Immune_Characteristics_sub <- Immune_Characteristics[Immune_Characteristics$TCGA.Study == i,]
  saveRDS(Immune_Characteristics_sub, file = paste0(i,"_The_Immune_Landscape_of_Cancer_Immune_Characteristics.rds"))
}






















# 2.test.R
COAD <- readRDS("COAD_The_Immune_Landscape_of_Cancer_Immune_Characteristics.rds")
colnames(COAD)
rowSums(COAD[,59:64])
rowSums(COAD[,5:7])









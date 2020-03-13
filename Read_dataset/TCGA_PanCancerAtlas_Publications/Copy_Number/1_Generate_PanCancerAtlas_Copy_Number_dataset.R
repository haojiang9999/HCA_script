#### 1_Generate_PanCancerAtlas_Copy_Number_dataset.R
### 1.Read_PanCancerAtlas_Copy_Number
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/00a32f7a-c85f-4f86-850d-be53973cbc4d/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg"
PANCAN_Genome_Wide_SNP_6<- read.table(file = filePath, sep = '\t', header = TRUE)
PANCAN_Genome_Wide_SNP_6[1:5,1:5]

























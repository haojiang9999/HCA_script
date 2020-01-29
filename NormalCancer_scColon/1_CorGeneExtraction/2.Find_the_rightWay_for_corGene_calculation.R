#### Find correlation details
source("/data8t_4/JH/MyJobs/1_R_script/RCA/RCA_seperate_fuction.R")
## I used correlation methods
CorRes <- corMatrix(log2.Tang.colon.cancer.FPKM.500,Ref.Tang.Adult.colon)
#### The that I calculate the cor for each gene
# One cell to reference table 
oneCell <- log2.Tang.colon.cancer.FPKM.500[,1,drop = F]
head(oneCell)
source("/data8t_4/JH/MyJobs/1_R_script/RCA/Step1_CorValueByGenes_1cell_vs_multiple_reference_cell_types.R")
CorRes.cell <- FindCorValueByGenes(oneCell,Ref.Tang.Adult.colon)
allGenes <- CorRes.cell$allGenes
BothAbove0 <- CorRes.cell$BothAbove0
singleAbove0 <- CorRes.cell$singleAbove0

## Sum of each gene cor
allGenesCor <- unlist(lapply(allGenes, sum)) ### This one was fit with corMatrix
BothAbove0Cor <- unlist(lapply(BothAbove0, sum))
singleAbove0Cor <- unlist(lapply(singleAbove0, sum))
CorCampare <- data.frame(allGenesCor = allGenesCor,
           BothAbove0Cor = BothAbove0Cor,
           singleAbove0Cor = singleAbove0Cor)


#### Conclusion  ######
# allGenesCor was the right corGene for my need













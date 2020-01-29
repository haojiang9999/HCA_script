#### Extract cor for each genes
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/NormalCancer_Core_script.R")
### Input was expression table and reference table
CorGene.res <- CorValueByGenes(log2.Tang.colon.cancer.FPKM.500,Ref.Tang.Adult.colon)
#test <- CorGene.res[["GSM2696954_scTrioSeq2Rna_CRC01_LN1_161"]]
#colSums(test)

CorValueResaults = CorGene.res
CorGene.Top.res <- TopCorGenes(CorGene.res)
CorGene.Top.res[["GSM2696954_scTrioSeq2Rna_CRC01_LN1_161"]]

#### 7.Immune cell expression
COAD.RSEM.gene.tpm
summary(colSums(COAD.RSEM.gene.tpm))
rownames(COAD.RSEM.gene.tpm)
# ENSG00000120217 = CD274/PDL1  ENSG00000197646= 	PDCD1LG2/PDL2
PDL1 <- COAD.RSEM.gene.tpm[grep("ENSG00000120217",rownames(COAD.RSEM.gene.tpm)),]
PDL2 <- COAD.RSEM.gene.tpm[grep("ENSG00000197646",rownames(COAD.RSEM.gene.tpm)),]
TCGA.COAD.PDLexp <- data.frame(PDL1.tpm = PDL1, PDL2.tpm = PDL2)
TCGA.COAD.PDLexp$rownames <- names(PDL1)






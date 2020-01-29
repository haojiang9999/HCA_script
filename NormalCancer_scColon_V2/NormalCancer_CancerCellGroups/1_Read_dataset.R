#### Read data
#1.Cancer cell groups
Cancer.group.list <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_scColon_V2/Cancer.group.list.exp.rds")
#2.TCGA COAD and TPM had been log2(x+0.001) transformed
#### read TCGA_COAD data ####
COAD_UCSC_Toil_tpm_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset.rds")
COAD.RSEM.gene.tpm <- COAD_UCSC_Toil_tpm_dataset$COAD.RSEM.gene.tpm
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset$COAD.pheno
gencode.v23.annotation <- COAD_UCSC_Toil_tpm_dataset$gencode.v23.annotation
head(gencode.v23.annotation)
## convert ensembleID to symbol
geneMatch <- match(rownames(COAD_RSEM_gene_tpm),gencode.v23.annotation$V1)
geneSymbol <- as.character(gencode.v23.annotation[geneMatch,]$V2)
COAD_tpm_symbol <- COAD_RSEM_gene_tpm
rownames(COAD_tpm_symbol) <- geneSymbol
head(COAD_tpm_symbol)
#3.Read sc Normal reference panel
scReference.V1 <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step2_Merge_and_Filter_the_scReference/2019_11_16_scReference.V1.rds")
scReference.list.CV.2000 <- scReference.V1$scReference.list.CV.2000









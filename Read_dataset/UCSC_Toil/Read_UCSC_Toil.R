#### 1.read RNA-seq expression data sets ####
#https://xenabrowser.net/datapages/?dataset=tcga_RSEM_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/gene_expression_RNAseq/tcga_RSEM_gene_tpm/tcga_RSEM_gene_tpm.gz"
#### Data (file names: *.rsem_genes.results) are downloaded, 
#    tpm values are extracted, log2(x+0.001) transformed
tcga_RSEM_gene_tpm <- read.table(filePath, header = T)
library(jsonlite)
tcga_RSEM_gene_tpm.metadata <- read_json("/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/gene_expression_RNAseq/tcga_RSEM_gene_tpm/tcga_RSEM_gene_tpm.json")
## add first column to rownames
rownames(tcga_RSEM_gene_tpm) <- tcga_RSEM_gene_tpm$sample
# remove first columns
tcga_RSEM_gene_tpm <- tcga_RSEM_gene_tpm[,-1]
head(tcga_RSEM_gene_tpm[1:10,1:10])
saveRDS(tcga_RSEM_gene_tpm, file = "tcga_RSEM_gene_tpm.rds")
#### 2.read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
COAD_Index <- phenotype_sp$cancer.type.abbreviation == "COAD"
COAD.pheno <- phenotype_sp[COAD_Index,]
head(COAD.pheno)
#### 3.Subset cancer types#### 
COAD.sampleID <- as.character(COAD.pheno$sample)
COAD.sampleID <- gsub("-",".",COAD.sampleID)
rownames(COAD.pheno) <- COAD.sampleID
## COAD sampleID in expression dataframe
COAD.sampleID.exp <- colnames(tcga_RSEM_gene_tpm) %in% COAD.sampleID
COAD_RSEM_gene_tpm  <- tcga_RSEM_gene_tpm[,COAD.sampleID.exp]
head()
COAD.pheno.exp <- COAD.pheno[colnames(COAD_RSEM_gene_tpm),]
## reverse log2 + 0.001 expression
COAD_RSEM_gene_tpm<-apply(COAD_RSEM_gene_tpm,2, function(x){
                          2^x - 0.001
                          })
# replace negative value to zeros
COAD_RSEM_gene_tpm[COAD_RSEM_gene_tpm <0 ] <- 0 
summary(colSums(COAD_RSEM_gene_tpm))

#### 4.gencode.v23.annotation.gene.probemap
geneAnnoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/gene_expression_RNAseq/tcga_RSEM_gene_tpm/gencode.v23.annotation.gene.probemap"
gencode.v23.annotation <- read.table(geneAnnoFilePath)
## biuld COAD data sets
COAD_UCSC_Toil_tpm <- list(COAD_RSEM_gene_tpm = COAD_RSEM_gene_tpm,
                           COAD.pheno.exp = COAD.pheno.exp,
                           gencode.v23.annotation = gencode.v23.annotation)
saveRDS(COAD_UCSC_Toil_tpm, file = "COAD_UCSC_Toil_tpm_dataset.rds")

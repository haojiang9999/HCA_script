#### read data ####
# Xena
#### 1.Read TCGA 450K methylation data sets ####
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/DNA_methylation/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz"
PANCAN_HumanMethylation450 <- read.delim(filePath, header = T)
saveRDS(PANCAN_HumanMethylation450, "PANCAN_HumanMethylation450.rds")
## add first column to rownames
rownames(PANCAN_HumanMethylation450) <- PANCAN_HumanMethylation450$sample
# remove first columns
PANCAN_HumanMethylation450 <- PANCAN_HumanMethylation450[,-1]
head(PANCAN_HumanMethylation450[1:10,1:10])
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
#### 3.Methylation annotation files
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/DNA_methylation/illuminaMethyl450_hg19_GPL16304_TCGAlegacy"
illuminaMethyl450_hg19_GPL16304_TCGAlegacy <- read.delim(filePath)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))

#### 4.Generate TCGA_Methyl450_hg19__dataset ####
#i="BRCA"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- phenotype_sp$cancer.type.abbreviation == i
  TCGA.pheno <- phenotype_sp[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$sample)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  rownames(TCGA.pheno) <- TCGA.sampleID
  ## Step3 Extract expression data by pheno data
  TCGA.sampleID.meth <- colnames(PANCAN_HumanMethylation450) %in% TCGA.sampleID
  TCGA_Methylation450  <- PANCAN_HumanMethylation450[,TCGA.sampleID.meth, drop = F]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_Methylation450),]
  ## Step4 Reverse log2 + 0.001 expression
  #TCGA_RSEM_gene_tpm<-apply(TCGA_RSEM_gene_tpm,2, function(x){
  #  2^x - 0.001
  #})
  # replace negative value to zeros
  #TCGA_RSEM_gene_tpm[TCGA_RSEM_gene_tpm <0 ] <- 0 
  ## Step5 Biuld TCGA data sets
  TCGA_Methylation450_dataset <- list(Methylation450 = TCGA_Methylation450,
                             pheno = TCGA.pheno.exp,
                             Methyl450_hg19_GPL16304 = illuminaMethyl450_hg19_GPL16304_TCGAlegacy)
  names(TCGA_Methylation450_dataset)<-c(paste0(i,"Methylation450"),pheno = paste0(i,".pheno"), "Methyl450_hg19_GPL16304")
  saveRDS(TCGA_Methylation450_dataset, file = paste0(i,"_PANCAN_Methylation450_dataset.rds"))
}
#TCGA_Methylation450[1:5,1:5]

#PANCAN_HumanMethylation450[1:5,1:5]

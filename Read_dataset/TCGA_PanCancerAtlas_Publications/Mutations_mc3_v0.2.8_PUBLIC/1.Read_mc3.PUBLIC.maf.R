#### 1.Read_mc3.PUBLIC.maf.R
### Step1.read data using maftools
library(maftools)
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/1c8cfe5f-e52d-41ba-94da-f15ea1337efc/mc3.v0.2.8.PUBLIC.maf"
mc3.PUBLIC.maf <- read.maf(maf = filePath)
saveRDS(mc3.PUBLIC.maf, file = "mc3.v0.2.8.PUBLIC.maf.rds")
mc3.v0.2.8.PUBLIC.maf <- readRDS("mc3.v0.2.8.PUBLIC.maf.rds")
## Metadata
mc3.v0.2.8.PUBLIC.metadata <- list(paper = "Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines",
                                   project = "PanCanAtlas Publications",
                                   url="https://gdc.cancer.gov/about-data/publications/pancanatlas")
#Shows sample summry.
getSampleSummary(mc3.v0.2.8.PUBLIC.maf)
#Shows all fields in MAF
getFields(mc3.v0.2.8.PUBLIC.maf)
# Gene mutation 
oncostrip(maf = mc3.v0.2.8.PUBLIC.maf, genes = c('KRAS','BRAF', 'APC','TP53'))
oncoplot(maf = mc3.v0.2.8.PUBLIC.maf, top = 10)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = mc3.v0.2.8.PUBLIC.maf, basename = 'mc3.v0.2.8.PUBLIC')
### Step2.Read summary files
TCGA.mc3.PUBLIC_geneSummary <- read.table("mc3.v0.2.8.PUBLIC_geneSummary.txt",header = T)
TCGA.mc3.PUBLIC_sampleSummary <- read.table("mc3.v0.2.8.PUBLIC_sampleSummary.txt",header = T)
TCGA.mc3.PUBLIC_summary <- read.table("mc3.v0.2.8.PUBLIC_summary.txt",header = T)
## change TCGA.mc3.PUBLIC_sampleSummary rownames
sampleID <- as.character(TCGA.mc3.PUBLIC_sampleSummary$Tumor_Sample_Barcode)
sampleID <- substr(sampleID,1,nchar(sampleID)-13)
rownames(TCGA.mc3.PUBLIC_sampleSummary) <- gsub("-",".",sampleID)
### Step3.Mutation genes by samples
# whether to include sysnonymous variants in ouput matrix. Default FALSE
TCGA.mc3.PUBLIC_mutCountMatrix <- mutCountMatrix(mc3.v0.2.8.PUBLIC.maf, removeNonMutated = F,
                                                 includeSyn = FALSE)
## Change sample names
sampleID <- colnames(TCGA.mc3.PUBLIC_mutCountMatrix)
sampleID <- substr(sampleID,1,nchar(sampleID)-13)
colnames(TCGA.mc3.PUBLIC_mutCountMatrix) <- gsub("-",".",sampleID)
TCGA.mc3.PUBLIC_mutCountMatrix[1:5,1:5]
### Step4.Phenotype data
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### Step5.Generate TCGA_UCSC_MC3_dataset ####
#i="ACC"

for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- phenotype_sp$cancer.type.abbreviation == i
  TCGA.pheno <- phenotype_sp[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$sample)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  rownames(TCGA.pheno) <- TCGA.sampleID
  ## Step3 Extract expression data by pheno data
  ####### TCGA.mc3.PUBLIC_sampleSummary ######
  TCGA.sampleID.exp <- rownames(TCGA.mc3.PUBLIC_sampleSummary) %in% TCGA.sampleID
  TCGA_tb_sampleSummary  <- TCGA.mc3.PUBLIC_sampleSummary[TCGA.sampleID.exp,]
  TCGA.pheno.exp <- TCGA.pheno[rownames(TCGA_tb_sampleSummary),]
  ####### TCGA.mc3.PUBLIC_mutCountMatrix ######
  TCGA.sampleID.exp <- colnames(TCGA.mc3.PUBLIC_mutCountMatrix) %in% TCGA.sampleID
  TCGA_tb_mutCountMatrix  <- TCGA.mc3.PUBLIC_mutCountMatrix[,TCGA.sampleID.exp]
  #TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA.mc3.PUBLIC_mutCountMatrix),]
  ## Step5 Biuld TCGA data sets
  TCGA_UCSC_dataset <- list(TCGA_tb_sampleSummary = TCGA_tb_sampleSummary,
                             TCGA_tb_mutCountMatrix =TCGA_tb_mutCountMatrix,
                             pheno = TCGA.pheno.exp,
                             mc3.v0.2.8.PUBLIC.metadata = mc3.v0.2.8.PUBLIC.metadata)
  names(TCGA_UCSC_dataset)<-c(paste0(i,".mc3.PUBLIC_sampleSummary"),paste0(i,".mc3.PUBLIC_mutCountMatrix"),pheno = paste0(i,".pheno"), "mc3.v0.2.8.PUBLIC.metadata")
  saveRDS(TCGA_UCSC_dataset, file = paste0(i,"_mc3.v0.2.8.PUBLIC.maf_dataset.rds"))
}















#### 1_Generate_Xena_TCGA_pancancer_Gene_level_non_silent_mutation_MC3.R
## 
#### 1.read mc3.v0.2.8.PUBLIC.nonsilentGene.xena ####
## TCGA Unified Ensemble "MC3" gene-level mutation calls. 1: non-silent mutation 0: wt
## https://xenabrowser.net/datapages/?dataset=mc3.v0.2.8.PUBLIC.nonsilentGene.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/somatic_mutation/Gene_level_non_silent_mutation/"
mc3.v0.2.8.PUBLIC.nonsilentGene.xena <- read.delim(paste0(filePath,"mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz"),
                                                                header = T)
mc3.v0.2.8.PUBLIC.nonsilentGene.xena[1:5,1:5]
class(mc3.v0.2.8.PUBLIC.nonsilentGene.xena)
mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx <- mc3.v0.2.8.PUBLIC.nonsilentGene.xena
GeneNames <- as.character(mc3.v0.2.8.PUBLIC.nonsilentGene.xena$sample)
## There was one NA name!!! why?
GeneNames[is.na(GeneNames)] <- "No_Name"
rownames(mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx) <- GeneNames
mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx <- mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx[,-1]
mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx[1:5,1:5]
## 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
## 3.Read probe Information
hugo_gencode_good_hg19_V24lift37_probemap <- read.table(file=paste0(filePath,"hugo_gencode_good_hg19_V24lift37_probemap"),header = T)
hugo_gencode_good_hg19_V24lift37_probemap[1:5,1:5]
## 4. Read metadata ####
library(rjson)
mc3.nonsilentGene.metadata <- fromJSON(file=paste0(filePath,"mc3.v0.2.8.PUBLIC.nonsilentGene.xena.json"))
## 5.Generate_Xena_TCGA_pancancer_PanCan33_Gene_level_non_silent_mutation_MC3_dataset ####
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
  TCGA.sampleID.exp <- colnames(mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx) %in% TCGA.sampleID
  TCGA_tb  <- mc3.v0.2.8.PUBLIC.nonsilentGene.xena_mx[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_tb),]
  ## Step4 Biuld TCGA data sets
  TCGA_PanCan33_mc3.nonsilentGene <- list(mc3.nonsilentGene = TCGA_tb,
                                    pheno = TCGA.pheno.exp,
                                    TCGA_tb.metadata = mc3.nonsilentGene.metadata,
                                    hugo_gencode_good_hg19_V24lift37_probemap = hugo_gencode_good_hg19_V24lift37_probemap)
  names(TCGA_PanCan33_mc3.nonsilentGene)<-c(paste0(i,".mc3.nonsilentGene.xena"),pheno = paste0(i,".pheno"), "mc3.nonsilentGene.metadata","Gene.Anno")
  saveRDS(TCGA_PanCan33_mc3.nonsilentGene, file = paste0(i,"_TCGA_PanCan33_mc3_nonsilentGene_dataset.rds"))
}



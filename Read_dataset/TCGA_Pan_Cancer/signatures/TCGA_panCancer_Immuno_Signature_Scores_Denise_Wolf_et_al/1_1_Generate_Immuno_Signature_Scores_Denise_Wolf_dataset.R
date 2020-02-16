#### 1_1_Generate_Immuno_Signature_Scores_Denise_Wolf_dataset.R
#### 1.read TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz ####
#### https://xenabrowser.net/datapages/?dataset=TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/signatures/Immuno_Signature_Scores_Denise_Wolf_et_al/"
TCGA_pancancer_Denise_Wolf_68ImmuneSigs <- read.table(paste0(filePath,"TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz"),
                                                      header = T, row.names = 1)
TCGA_pancancer_Denise_Wolf_68ImmuneSigs[1:5,1:5]
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read metadata ####
library(rjson)
Denise_Wolf_68ImmuneSigs.metadata <- fromJSON(file=paste0(filePath,"TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.json"))
#### 4.Generate_Xena_TCGA_pancancer_Subtype_Immune_dataset ####
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
  TCGA.sampleID.exp <- colnames(TCGA_pancancer_Denise_Wolf_68ImmuneSigs) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- TCGA_pancancer_Denise_Wolf_68ImmuneSigs[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Biuld TCGA data sets
  TCGA_Denise_Wolf_68ImmuneSigs.xena <- list(TCGA_pancancer_Denise_Wolf_68ImmuneSigs = TCGA_RSEM_gene,
                                   pheno = TCGA.pheno.exp,
                                   Denise_Wolf_68ImmuneSigs.metadata = Denise_Wolf_68ImmuneSigs.metadata)
  names(TCGA_Denise_Wolf_68ImmuneSigs.xena)<-c(paste0(i,".pancancer.Denise.Wolf.68ImmuneSigs.xena"),pheno = paste0(i,".pheno"), "Denise_Wolf_68ImmuneSigs.metadata")
  saveRDS(TCGA_Denise_Wolf_68ImmuneSigs.xena, file = paste0(i,"_Denise_Wolf_68ImmuneSigs_dataset.rds"))
}
























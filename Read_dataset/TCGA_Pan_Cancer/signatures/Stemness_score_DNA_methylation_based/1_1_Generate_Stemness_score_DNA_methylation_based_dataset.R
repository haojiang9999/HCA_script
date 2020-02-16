# 1_1_Generate_Stemness_score_DNA_methylation_based_dataset.R
#### 1.read StemnessScores_DNAmeth_20170210.tsv ####
#### https://xenabrowser.net/datapages/?dataset=StemnessScores_DNAmeth_20170210.tsv&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/signatures/Stemness_score_DNA_methylation_based/"
TCGA_pancancer_StemnessScores_DNAmeth <- read.table(paste0(filePath,"StemnessScores_DNAmeth_20170210.tsv.gz"),
                                                      header = T, row.names = 1)
# DNAss: DNA methylation-based (Stem cell signature probes (219 probes), that combines the 3 signatures listed below). This score will drive the main figures in the PancanAtlas paper
# EREG-METHss: Epigenetically regulated DNA methylation-based (87 probes).
# DMPss: Differentially methylated probes-based (62 probes)
# ENHss: Enhancer Elements/DNAmethylation-based (82 probes)
TCGA_pancancer_StemnessScores_DNAmeth[1:4,1:4]
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read metadata ####
library(rjson)
StemnessScores_DNAmeth.metadata <- fromJSON(file=paste0(filePath,"StemnessScores_DNAmeth_20170210.tsv.json"))
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
  TCGA.sampleID.exp <- colnames(TCGA_pancancer_StemnessScores_DNAmeth) %in% TCGA.sampleID
  TCGA_df  <- TCGA_pancancer_StemnessScores_DNAmeth[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_df),]
  ## Step4 Biuld TCGA data sets
  TCGA_StemnessScores_DNAmeth.xena <- list(TCGA_pancancer_StemnessScores_DNAmeth = TCGA_df,
                                             pheno = TCGA.pheno.exp,
                                             StemnessScores_DNAmeth.metadata = StemnessScores_DNAmeth.metadata)
  names(TCGA_StemnessScores_DNAmeth.xena)<-c(paste0(i,".pancancer.StemnessScores.DNAmeth.xena"),pheno = paste0(i,".pheno"), "StemnessScores_DNAmeth.metadata")
  saveRDS(TCGA_StemnessScores_DNAmeth.xena, file = paste0(i,"_StemnessScores_DNAmeth_dataset.rds"))
}



















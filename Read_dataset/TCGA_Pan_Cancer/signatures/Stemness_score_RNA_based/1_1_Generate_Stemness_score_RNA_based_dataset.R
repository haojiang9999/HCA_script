# 1_1_Generate_Stemness_score_RNA_based_dataset.R
#### 1.read StemnessScores_RNAexp_20170127.2.tsv.gz ####
#### https://xenabrowser.net/datapages/?dataset=StemnessScores_RNAexp_20170127.2.tsv&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/signatures/Stemness_score_RNA_based_n_10876/"
TCGA_pancancer_StemnessScores_RNAexp <- read.table(paste0(filePath,"StemnessScores_RNAexp_20170127.2.tsv.gz"),
                                                   header = T, row.names = 1)
# RNAss: RNA expression-based (All set of available genes) this score will drive the main figures in PancanAtlas paper.
# EREG.EXPss: Epigenetically regulated RNA expression-based (103 genes).
TCGA_pancancer_StemnessScores_RNAexp[1:2,1:5] 
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read metadata ####
library(rjson)
StemnessScores_RNAexp.metadata <- fromJSON(file=paste0(filePath,"StemnessScores_RNAexp_20170127.2.tsv.json"))
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
  TCGA.sampleID.exp <- colnames(TCGA_pancancer_StemnessScores_RNAexp) %in% TCGA.sampleID
  TCGA_df  <- TCGA_pancancer_StemnessScores_RNAexp[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_df),]
  ## Step4 Biuld TCGA data sets
  TCGA_StemnessScores_RNAexp.xena <- list(TCGA_pancancer_StemnessScores_RNAexp = TCGA_df,
                                           pheno = TCGA.pheno.exp,
                                          StemnessScores_RNAexp.metadata = StemnessScores_RNAexp.metadata)
  names(TCGA_StemnessScores_RNAexp.xena)<-c(paste0(i,".pancancer.StemnessScores.RNAexp.xena"),pheno = paste0(i,".pheno"), "StemnessScores_RNAexp.metadata")
  saveRDS(TCGA_StemnessScores_RNAexp.xena, file = paste0(i,"_StemnessScores_RNAexp_dataset.rds"))
}
                                                   
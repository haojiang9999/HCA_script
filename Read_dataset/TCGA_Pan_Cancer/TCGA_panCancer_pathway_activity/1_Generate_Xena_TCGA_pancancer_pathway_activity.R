#### 1_Generate_Xena_TCGA_pancancer_pathway_activity.R
##################ssGSEA data###############################
#### 1.read PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt ####
#### https://xenabrowser.net/datapages/?dataset=PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/pathway_activity/"
PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level <- read.delim(paste0(filePath,"PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt.gz"),
                                                                header = T)
PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level[1:5,1:5]
rownames(PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level) <- PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level$X
saveRDS(PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level, file = "PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.rds")
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read metadata ####
library(rjson)
PanCan33_ssGSEA.metadata <- fromJSON(file=paste0(filePath,"PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt.json"))
#### 4.Generate_Xena_TCGA_pancancer_PanCan33_ssGSEA_1387GeneSets_dataset ####
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
  TCGA.sampleID.exp <- colnames(PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Biuld TCGA data sets
  TCGA_PanCan33_ssGSEA.xena <- list(TCGA_PanCan33_ssGSEA.xena = TCGA_RSEM_gene,
                               pheno = TCGA.pheno.exp,
                               PanCan33_ssGSEA.metadata = PanCan33_ssGSEA.metadata)
  names(TCGA_PanCan33_ssGSEA.xena)<-c(paste0(i,".PanCan33.ssGSEA.xena"),pheno = paste0(i,".pheno"), "PanCan33_ssGSEA.metadata")
  saveRDS(TCGA_PanCan33_ssGSEA.xena, file = paste0(i,"_TCGA_PanCan33_ssGSEA_dataset.rds"))
}
























#### 1_1_Generate_Xena_TCGA_pancancer_Immune_subtype.R
#### 1.read Subtype_Immune_Model_Based.txt.gz ####
#### https://xenabrowser.net/datapages/?dataset=Subtype_Immune_Model_Based.txt&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
## Paper:Co-expression modules identified from published immune signatures reveal five distinct immune subtypes in breast cancer.
## doi: 10.1007/s10549-016-4041-3.
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/Phenotype/"
Subtype_Immune_Model_Based <- read.delim(paste0(filePath,"Subtype_Immune_Model_Based.txt.gz"))
TCGA.sampleID <- as.character(Subtype_Immune_Model_Based$sample)
TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
rownames(Subtype_Immune_Model_Based) <- TCGA.sampleID
Subtype_Immune_Model_Based$rownames <- TCGA.sampleID
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read metadata ####
library(rjson)
Subtype_Immune.metadata <- fromJSON(file=paste0(filePath,"Subtype_Immune_Model_Based.txt.json"))
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
  TCGA.sampleID.exp <- rownames(Subtype_Immune_Model_Based) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- Subtype_Immune_Model_Based[TCGA.sampleID.exp,]
  TCGA.pheno.exp <- TCGA.pheno[rownames(TCGA_RSEM_gene),]
  ## Step4 Biuld TCGA data sets
  TCGA_Subtype_Immune.xena <- list(TCGA_Subtype_Immune.xena = TCGA_RSEM_gene,
                                    pheno = TCGA.pheno.exp,
                                   Subtype_Immune.metadata = Subtype_Immune.metadata)
  names(TCGA_Subtype_Immune.xena)<-c(paste0(i,".Subtype.Immune.xena"),pheno = paste0(i,".pheno"), "Subtype_Immune.metadata")
  saveRDS(TCGA_Subtype_Immune.xena, file = paste0(i,"_Subtype_Immune_dataset.rds"))
}

                                         
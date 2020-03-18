#### 1.Generate_TCGA_Xena_survaival_dataset.R
# Paper:"An Integrated TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR) to drive high quality survival outcome analytics"
### 1.Read Survival_SupplementalTable_S1_20171025_xena_sp.gz
SurvivalFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
Survival_SupplementalTable <- read.delim(SurvivalFilePath) 
table(Survival_SupplementalTable$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(Survival_SupplementalTable$cancer.type.abbreviation))
### 2.Metadata
library(jsonlite)
Survival_SupplementalTable.metadata <- read_json("/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.json")

### 3.Generate_TCGA_Xena_survaival_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- Survival_SupplementalTable$cancer.type.abbreviation == i
  TCGA.pheno <- Survival_SupplementalTable[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$sample)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  rownames(TCGA.pheno) <- TCGA.sampleID
  TCGA.pheno$rownames <- TCGA.sampleID
  ## Step3 Biuld TCGA data sets
  Survival_SupplementalTable.xena <- list(TCGA.pheno = TCGA.pheno,Survival_SupplementalTable.metadata = Survival_SupplementalTable.metadata)
  names(Survival_SupplementalTable.xena)<-c(paste0(i,".Survival_SupplementalTable.xena"),pheno = paste0(i,".Survival_SupplementalTable.metadata"))
  saveRDS(Survival_SupplementalTable.xena, file = paste0(i,"_Survival_SupplementalTabl_dataset.rds"))
}



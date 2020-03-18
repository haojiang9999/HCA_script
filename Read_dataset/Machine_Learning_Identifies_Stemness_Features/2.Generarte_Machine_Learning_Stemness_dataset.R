#### 2.Generarte_Machine_Learning_Stemness_dataset.R
### 1.Metadata
Machine_Learning_Stemness.metadata <- list(paper = "Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation",
                                           ExtraData_url = "https://gdc.cancer.gov/about-data/publications/PanCanStemness-2018",
                                           immuneCell.Absolute.df = "For TotalLeukocyte was  ESTIMATE applied to DNA methylation data We also obtained absolute estimates by scaling their relative abundance by overall leukocyte infiltration in each tumor as determined by ESTIMATE applied to DNA methylation data 
                                           1)Extract immune cell info:But this is relative infiltration of each cell type
                                           2)Get the absolute infiltration value for immune cells:Relative value multiple TotalLeukocyte(ESTIMATE)")

### 2.Therre were more data down load from 
# https://gdc.cancer.gov/about-data/publications/PanCanStemness-2018
#TCGA-CDR-SupplementalTableS1 <- re("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/TCGA-CDR-SupplementalTableS1.xlsx")
TCGA.Kallisto.cibersort.relative <- read.delim("/stor/jianghao/Paper_data/Machine_Learning_Identifies_Stemness_Features/TCGA.Kallisto.cibersort.relative.tsv")
table(TCGA.Kallisto.cibersort.relative$CancerType)
TCGA.cancer.types <- names(table(TCGA.Kallisto.cibersort.relative$CancerType))

#### 4.Generate TCGA_UCSC_Toil_tpm_dataset ####
#i="ACC"

for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA.Kallisto.cibersort.relative$CancerType == i
  TCGA.Kallisto.cibersort.relative_sub <- TCGA.Kallisto.cibersort.relative[TCGA_Index,]
  #### Step2 convert SampleID to rownames
  #TCGA.Kallisto.cibersort.relative_sub$SampleID
  TCGA.Kallisto.cibersort.relative_sub$rownames <- substr(TCGA.Kallisto.cibersort.relative_sub$SampleID, start = 1, stop = 15)
  ### Step3 Convert to absolute version 
  #For "TotalLeukocyte" was  ESTIMATE applied to DNA methylation data
  #We also obtained absolute estimates by scaling their relative abundance by overall
  #leukocyte infiltration in each tumor as determined by ESTIMATE applied to DNA methylation data 
  ## 1)Extract immune cell info
  #colnames(TCGA.Kallisto.cibersort.relative_sub)
  immuneCell.df <- TCGA.Kallisto.cibersort.relative_sub[,3:24] 
  ## But this is relative infiltration of each cell type
  ## 2)Get the absolute infiltration value for immune cells
  # Relative value multiple TotalLeukocyte(ESTIMATE)
  immuneCell.Absolute.df <- as.data.frame(apply(immuneCell.df,2,function(x){x*TCGA.Kallisto.cibersort.relative_sub$TotalLeukocyte}))
  immuneCell.Absolute.df$rownames <- TCGA.Kallisto.cibersort.relative_sub$rownames
  
  ### Step4 Biuld TCGA data sets
  Machine_Learning_Stemness_dataset <- list(TCGA.Kallisto.cibersort.relative_sub = TCGA.Kallisto.cibersort.relative_sub,
                                            immuneCell.Absolute.df = immuneCell.Absolute.df,
                                            Machine_Learning_Stemness.metadata = Machine_Learning_Stemness.metadata)
  names(Machine_Learning_Stemness_dataset)<-c(paste0(i,"_Machine_Learning_Stemness.ImmuneCell.cibersort.relative"), 
                                              paste0(i,"_Machine_Learning_StemnessimmuneCell.Absolute"),
                                              "Machine_Learning_Stemness.metadata")
  saveRDS(Machine_Learning_Stemness_dataset, file = paste0(i,"_Machine_Learning_Stemness.ImmuneCell.cibersort.relative_dataset.rds"))
}

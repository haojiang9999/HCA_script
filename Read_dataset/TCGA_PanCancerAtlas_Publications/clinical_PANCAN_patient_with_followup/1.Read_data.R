#### 1.Read_data.R
clinical_PANCAN_patient_with_followup <- read.delim(file = 'clinical_PANCAN_patient_with_followup.tsv', sep = '\t')
clinical_PANCAN_patient_with_followup[1:5,1:5]
## TSS Tissue source site
TSS.info <- clinical_PANCAN_patient_with_followup[,2:3]
TSS.info$barcode <- gsub("-",".",TSS.info$bcr_patient_barcode)
### 2.Save data
## PanCanAtlas Publications
TSS.info.list <- list(TSS.info = TSS.info, meatadata = "https://gdc.cancer.gov/about-data/publications/pancanatlas")
saveRDS(TSS.info.list, file = "TCGA_PanCancerAtlas_public_TSS.info.rds")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
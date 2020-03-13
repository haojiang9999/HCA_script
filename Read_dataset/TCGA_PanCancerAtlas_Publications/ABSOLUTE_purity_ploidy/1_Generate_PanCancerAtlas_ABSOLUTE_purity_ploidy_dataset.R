#### 1_Generate_PanCancerAtlas_ABSOLUTE_purity_ploidy_dataset.R
### 1.Read_PanCancerAtlas_ABSOLUTE_purity_ploidy
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/4f277128-f793-4354-a13d-30cc7fe9f6b5/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
mastercalls.abs_tables<- read.table(file = filePath, sep = '\t', header = TRUE)
purity.ploidy <- mastercalls.abs_tables
purity.ploidy$rownames <- gsub("-",".",purity.ploidy$array)
rownames(purity.ploidy) <- gsub("-",".",purity.ploidy$array)
purity.ploidy[1:5,]
### 2.Pheno data
TCGA_PanCancerAtlas_CDR_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/clinical_PANCAN_patient_with_followup/TCGA_PanCancerAtlas_CDR_dataset.rds")
TCGA_CDR <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR
TCGA_CDR_Notes <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR_Notes
metadata <- TCGA_PanCancerAtlas_CDR_dataset$metadata
TCGA.cancer.types <- names(table(TCGA_CDR$type))

### 3.RNA sample barcode
sample.barcode <- rownames(purity.ploidy)
sample.barcode <- substr(sample.barcode, start = 1, stop = 12)
#### 4.Generate PanCancerAtlas_RNA_final_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA_CDR$type == i
  TCGA.pheno <- TCGA_CDR[TCGA_Index,]
  ## Step2 Extract expression data by pheno data
  TCGA.sampleID.exp <- rownames(purity.ploidy)[sample.barcode %in% TCGA.pheno$barcode]
  TCGA_tb  <- purity.ploidy[TCGA.sampleID.exp,]
  TCGA_tb_sampleID<- substr(rownames(TCGA_tb), start = 1, stop = 12)
  TCGA.pheno.exp <- TCGA.pheno[TCGA.pheno$barcode %in% TCGA_tb_sampleID,]
  ## Step3 Biuld TCGA data sets
  TCGA_PanCancerAtlas <- list(purity.ploidy = TCGA_tb,
                              pheno = TCGA.pheno.exp,
                              phenoNotes = TCGA_CDR_Notes,
                              metadata = metadata)
  names(TCGA_PanCancerAtlas)<-c(paste0(i,".purity.ploidy"),pheno = paste0(i,".pheno"), "phenoNotes","metadata")
  saveRDS(TCGA_PanCancerAtlas, file = paste0(i,"_PanCancerAtlas_Publish_purity_ploidy_dataset.rds"))
}












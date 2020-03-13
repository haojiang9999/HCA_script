#### 1_Generate_PanCancerAtlas_RPPA_final_dataset.R
### 1.Read RPPA data
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/fcbb373e-28d4-4818-92f3-601ede3da5e1/TCGA-RPPA-pancan-clean.txt"
RPPA_final<- read.table(file = filePath, sep = '\t', header = TRUE)
RPPA_final[1:5,1:5]
SampleID <- RPPA_final$SampleID
SampleID <- gsub("-",".",SampleID)
rownames(RPPA_final) <- substr(SampleID, start = 1, stop = 15)
RPPA_final_mx <- RPPA_final
RPPA_final_mx[1:5,1:5]
RPPA_final_mx <- RPPA_final_mx[,-c(1:2)]
RPPA_final_mx <- as.matrix(RPPA_final_mx)
RPPA_final_mx <- t(RPPA_final_mx)
### 2.Pheno data
TCGA_PanCancerAtlas_CDR_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/clinical_PANCAN_patient_with_followup/TCGA_PanCancerAtlas_CDR_dataset.rds")
TCGA_CDR <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR
TCGA_CDR_Notes <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR_Notes
metadata <- TCGA_PanCancerAtlas_CDR_dataset$metadata
TCGA.cancer.types <- names(table(TCGA_CDR$type))
table(RPPA_final$TumorType)
#### 4.Generate PanCancerAtlas_RNA_final_dataset ####
#i="COAD"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- rownames(RPPA_final[RPPA_final$TumorType == i,])
  ## Step2 Extract expression data by pheno data
  TCGA_tb  <- RPPA_final_mx[,TCGA_Index]
  TCGA_tb_sampleID<- substr(colnames(TCGA_tb), start = 1, stop = 12)
  TCGA.pheno.exp <- TCGA_CDR[TCGA_CDR$barcode %in% TCGA_tb_sampleID,]
  ## Step3 Biuld TCGA data sets
  TCGA_PanCancerAtlas <- list(RPPA_final = TCGA_tb,
                              pheno = TCGA.pheno.exp,
                              phenoNotes = TCGA_CDR_Notes,
                              metadata = metadata)
  names(TCGA_PanCancerAtlas)<-c(paste0(i,".RPPA.final.mx"),pheno = paste0(i,".pheno"), "phenoNotes","metadata")
  saveRDS(TCGA_PanCancerAtlas, file = paste0(i,"_PanCancerAtlas_Publish_RPPA_Final_dataset.rds"))
}






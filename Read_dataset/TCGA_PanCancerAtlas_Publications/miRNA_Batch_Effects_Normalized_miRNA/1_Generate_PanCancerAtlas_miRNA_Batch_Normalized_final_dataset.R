#### 1_Generate_PanCancerAtlas_miRNA_Batch_Normalized_final_dataset.R
### 1.Read_PanCancerAtlas_RNA_final
info.filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/55d9bf6f-0712-4315-b588-e6f8e295018e/PanCanAtlas_miRNA_sample_information_list.txt"
PanCanAtlas_miRNA_sample_information_list<- read.table(file = info.filePath, sep = '\t', header = TRUE)
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/1c6174d9-8ffb-466e-b5ee-07b204c15cf8/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv"
pancanMiRs <-  read.csv(file = filePath, header = TRUE)
pancanMiRs[1:5,1:5]
Genes <- pancanMiRs$Genes
pancanMiRs_mx <- pancanMiRs[,-c(1,2)]
pancanMiRs_mx[1:5,1:5]
rownames(pancanMiRs_mx) <- Genes
sampleID <- colnames(pancanMiRs_mx)
colnames(pancanMiRs_mx)<- substr(sampleID, start = 1, stop = 15)
pancanMiRs_mx[1:5,1:5]
### 2.Pheno data
TCGA_PanCancerAtlas_CDR_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/clinical_PANCAN_patient_with_followup/TCGA_PanCancerAtlas_CDR_dataset.rds")
TCGA_CDR <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR
TCGA_CDR_Notes <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR_Notes
metadata <- TCGA_PanCancerAtlas_CDR_dataset$metadata
TCGA.cancer.types <- names(table(TCGA_CDR$type))
### 3.RNA sample barcode
sample.barcode <- colnames(pancanMiRs_mx)
sample.barcode <- substr(sample.barcode, start = 1, stop = 12)

#### 4.Generate PanCancerAtlas_RNA_final_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA_CDR$type == i
  TCGA.pheno <- TCGA_CDR[TCGA_Index,]
  ## Step2 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(pancanMiRs_mx)[sample.barcode %in% TCGA.pheno$barcode]
  TCGA_tb  <- pancanMiRs_mx[,TCGA.sampleID.exp]
  TCGA_tb_sampleID<- substr(colnames(TCGA_tb), start = 1, stop = 12)
  TCGA.pheno.exp <- TCGA.pheno[TCGA.pheno$barcode %in% TCGA_tb_sampleID,]
  ## Step3 Biuld TCGA data sets
  TCGA_PanCancerAtlas <- list(pancanMiRs_mx = TCGA_tb,
                              pheno = TCGA.pheno.exp,
                              phenoNotes = TCGA_CDR_Notes,
                              sample_information = PanCanAtlas_miRNA_sample_information_list,
                              metadata = metadata)
  names(TCGA_PanCancerAtlas)<-c(paste0(i,".pancan.miRNA.mx"),pheno = paste0(i,".pheno"), "phenoNotes","sample_information","metadata")
  saveRDS(TCGA_PanCancerAtlas, file = paste0(i,"_PanCancerAtlas_Publish_miRNA_Batch_Normalized_dataset.rds"))
}
























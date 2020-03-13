#### 1_Generate_PanCancerAtlas_PARADIGM_Pathway_dataset.R
### 1.Read_PanCancerAtlas_PARADIGM_Pathway
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/7d4c0344-f018-4ab0-949a-09815f483480/merge_merged_reals.tar.gz"
untar(filePath,files="merge_merged_reals.txt")
PARADIGM_Pathway_Inference_Matrix<- read.table("merge_merged_reals.txt", sep = '\t', header = TRUE)
PARADIGM_Pathway_Inference_Matrix[1:5,1:5]
pathway <- as.character(PARADIGM_Pathway_Inference_Matrix$Gene)
PARADIGM_mx <- PARADIGM_Pathway_Inference_Matrix[,-1]
PARADIGM_mx[1:5,1:5]
rownames(PARADIGM_mx) <- pathway
PARADIGM_mx <- as.matrix(PARADIGM_mx)
saveRDS(PARADIGM_mx, file = "All_PanCancerAtlas_Publish_PARADIGM_pathway_dataset.rds")
### 2.Pheno data
TCGA_PanCancerAtlas_CDR_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/clinical_PANCAN_patient_with_followup/TCGA_PanCancerAtlas_CDR_dataset.rds")
TCGA_CDR <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR
TCGA_CDR_Notes <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR_Notes
metadata <- TCGA_PanCancerAtlas_CDR_dataset$metadata
TCGA.cancer.types <- names(table(TCGA_CDR$type))

#### 4.Generate PanCancerAtlas_RNA_final_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA_CDR$type == i
  TCGA.pheno <- TCGA_CDR[TCGA_Index,]
  ## Step2 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(PARADIGM_mx)[colnames(PARADIGM_mx) %in% TCGA.pheno$barcode]
  TCGA_tb  <- PARADIGM_mx[,TCGA.sampleID.exp]
  TCGA_tb_sampleID<- substr(colnames(TCGA_tb), start = 1, stop = 12)
  TCGA.pheno.exp <- TCGA.pheno[TCGA.pheno$barcode %in% TCGA_tb_sampleID,]
  ## Step3 Biuld TCGA data sets
  TCGA_PanCancerAtlas <- list(PARADIGM_mx = TCGA_tb,
                              pheno = TCGA.pheno.exp,
                              phenoNotes = TCGA_CDR_Notes,
                              pathway = pathway,
                              metadata = metadata)
  names(TCGA_PanCancerAtlas)<-c(paste0(i,".PARADIGM.pathway.mx"),pheno = paste0(i,".pheno"), "phenoNotes","pathway","metadata")
  saveRDS(TCGA_PanCancerAtlas, file = paste0(i,"_PanCancerAtlas_Publish_PARADIGM_pathway_dataset.rds"))
}

#### 1_Generate_PanCancerAtlas_RNA_final_dataset.R
### 1.Read_PanCancerAtlas_RNA_final
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/3586c0da-64d0-4b74-a449-5ff4d9136611/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
RNA_final<- read.table(file = filePath, sep = '\t', header = TRUE)
saveRDS(RNA_final, file = "RNA_final.rds")
colnames(RNA_final)
table(is.na(RNA_final))
## It's real merge all gene expression like array and RNA-seq
dim(RNA_final)
RNA_final[1:5,1:5]
geneAnno <- as.character(RNA_final$gene_id)
geneAnno<- data.frame(GeneSymbol = unlist(lapply(strsplit(geneAnno,split='|', fixed=TRUE), '[[', 1)),
           GeneID = unlist(lapply(strsplit(geneAnno,split='|', fixed=TRUE), '[[', 2)))
RNA_final_mx <- RNA_final[,-1]
class(RNA_final)
RNA_final_mx[1:5,1:5]
RNA_final_mx <- as.matrix(RNA_final_mx)
rownames(RNA_final_mx) <- geneAnno$GeneSymbol
sampleID <- colnames(RNA_final_mx)
colnames(RNA_final_mx)<- substr(sampleID, start = 1, stop = 15)
RNA_final_mx[1:40,1:5]
table(is.na(RNA_final_mx))

### 2.Pheno data
TCGA_PanCancerAtlas_CDR_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/clinical_PANCAN_patient_with_followup/TCGA_PanCancerAtlas_CDR_dataset.rds")
TCGA_CDR <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR
TCGA_CDR_Notes <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR_Notes
metadata <- TCGA_PanCancerAtlas_CDR_dataset$metadata
TCGA.cancer.types <- names(table(TCGA_CDR$type))
### 3.RNA sample barcode
sample.barcode <- colnames(RNA_final_mx)
sample.barcode <- substr(sample.barcode, start = 1, stop = 12)
#### 4.Generate PanCancerAtlas_RNA_final_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA_CDR$type == i
  TCGA.pheno <- TCGA_CDR[TCGA_Index,]
  ## Step2 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(RNA_final_mx)[sample.barcode %in% TCGA.pheno$barcode]
  TCGA_tb  <- RNA_final_mx[,TCGA.sampleID.exp]
  TCGA_tb_sampleID<- substr(colnames(TCGA_tb), start = 1, stop = 12)
  TCGA.pheno.exp <- TCGA.pheno[TCGA.pheno$barcode %in% TCGA_tb_sampleID,]
  ## Step3 Biuld TCGA data sets
  TCGA_PanCancerAtlas <- list(RNA_final = TCGA_tb,
                             pheno = TCGA.pheno.exp,
                             phenoNotes = TCGA_CDR_Notes,
                             geneAnno = geneAnno,
                             metadata = metadata)
  names(TCGA_PanCancerAtlas)<-c(paste0(i,".RNA.final.mx"),pheno = paste0(i,".pheno"), "phenoNotes","geneAnno","metadata")
  saveRDS(TCGA_PanCancerAtlas, file = paste0(i,"_PanCancerAtlas_Publish_RNA_Final_dataset.rds"))
}


























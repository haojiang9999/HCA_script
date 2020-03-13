##### 1_Generate_PanCancerAtlas_merged_Methylation_dataset.R
### 1.Read_DNA Methylation (Merged 27K+450K Only)
filePath <- "/stor/jianghao/Paper_data/PanCancerAtlas_Publications/d82e2c44-89eb-43d9-b6d3-712732bf6a53/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv"
merged_HumanMethylation27_HumanMethylation450.betaValue<- read.table(file = filePath, sep = '\t', header = TRUE)
merged_HumanMethylation27_HumanMethylation450.betaValue[1:5,1:5]
mergedMethyl <- merged_HumanMethylation27_HumanMethylation450.betaValue[,-1]
cgID <- merged_HumanMethylation27_HumanMethylation450.betaValue$Composite.Element.REF
rownames(mergedMethyl) <- cgID
mergedMethyl[1:5,1:5]
colnames(mergedMethyl)<- substr(colnames(merged_HumanMethylation27_HumanMethylation450.betaValue[,-1]), start = 1, stop = 15)
mergedMethyl <- as.matrix(mergedMethyl)
class(mergedMethyl)

### 2.Pheno data
TCGA_PanCancerAtlas_CDR_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_PanCancerAtlas_Publications/clinical_PANCAN_patient_with_followup/TCGA_PanCancerAtlas_CDR_dataset.rds")
TCGA_CDR <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR
TCGA_CDR_Notes <- TCGA_PanCancerAtlas_CDR_dataset$TCGA_CDR_Notes
metadata <- TCGA_PanCancerAtlas_CDR_dataset$metadata
TCGA.cancer.types <- names(table(TCGA_CDR$type))
### 3.RNA sample barcode
sample.barcode <- colnames(mergedMethyl)
sample.barcode <- substr(sample.barcode, start = 1, stop = 12)

#### 4.Methylation annotation files
MethfilePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/DNA_methylation/illuminaMethyl450_hg19_GPL16304_TCGAlegacy"
illuminaMethyl450_hg19_GPL16304_TCGAlegacy <- read.delim(MethfilePath)
head(illuminaMethyl450_hg19_GPL16304_TCGAlegacy)
## all cgID been annotated
table(cgID %in% illuminaMethyl450_hg19_GPL16304_TCGAlegacy$X.id)
cgAnno <- illuminaMethyl450_hg19_GPL16304_TCGAlegacy[illuminaMethyl450_hg19_GPL16304_TCGAlegacy$X.id %in% cgID,]
#### 5.Generate PanCancerAtlas_RNA_final_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA_CDR$type == i
  TCGA.pheno <- TCGA_CDR[TCGA_Index,]
  ## Step2 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(mergedMethyl)[sample.barcode %in% TCGA.pheno$barcode]
  TCGA_tb  <- mergedMethyl[,TCGA.sampleID.exp]
  TCGA_tb_sampleID<- substr(colnames(TCGA_tb), start = 1, stop = 12)
  TCGA.pheno.exp <- TCGA.pheno[TCGA.pheno$barcode %in% TCGA_tb_sampleID,]
  ## Step3 Biuld TCGA data sets
  TCGA_PanCancerAtlas <- list(mergedMethyl = TCGA_tb,
                              pheno = TCGA.pheno.exp,
                              phenoNotes = TCGA_CDR_Notes,
                              cgAnno = cgAnno,
                              metadata = metadata)
  names(TCGA_PanCancerAtlas)<-c(paste0(i,".mergedMethyl.27.450.mx"),pheno = paste0(i,".pheno"), "phenoNotes","cgAnno","metadata")
  saveRDS(TCGA_PanCancerAtlas, file = paste0(i,"_PanCancerAtlas_Publish_mergedMethyl_27K_450K_dataset.rds"))
}








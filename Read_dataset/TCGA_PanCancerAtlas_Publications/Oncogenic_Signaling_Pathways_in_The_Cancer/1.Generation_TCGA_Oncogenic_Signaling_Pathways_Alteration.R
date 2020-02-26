#### 1.Generation_TCGA_Oncogenic_Signaling_Pathways_Alteration.R
# Paper:Oncogenic Signaling Pathways in The Cancer
### 1.Read Genomic Alteration Matrices_pathway_level table
Alteration_pathway_level <- read.csv("IM_Table S4. Genomic Alteration Matrices_pathway_level.csv")
sampleID <- as.character(Alteration_pathway_level$SAMPLE_BARCODE)
sampleID <- gsub("-",".",sampleID)
Alteration_pathway_level$rownames <- sampleID
rownames(Alteration_pathway_level)  <- sampleID
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Metadata generation ####
Alteration_pathway_level.metadata <- list(Paper = "Oncogenic Signaling Pathways in The Cancer",
                                          Table = "Table S4. Genomic Alteration Matrices_pathway_level.csv",
                                          doi =  "10.1016/j.cell.2018.03.035.")
#### 4.Generate_Oncogenic_Signaling_Pathways_Alteration_dataset ####
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
  TCGA.sampleID.exp <- Alteration_pathway_level$rownames %in% TCGA.sampleID
  TCGA_df  <- Alteration_pathway_level[TCGA.sampleID.exp,]
  ## Step4 Biuld TCGA data sets
  TCGA_Oncogenic_Signaling_Pathways_Alteration <- list(Alteration_pathway_level = TCGA_df,
                                          Alteration_pathway_level.metadata = Alteration_pathway_level.metadata)
  names(TCGA_Oncogenic_Signaling_Pathways_Alteration)<-c(paste0(i,".Alteration.pathway.level"), "Alteration_pathway_level.metadata")
  saveRDS(TCGA_Oncogenic_Signaling_Pathways_Alteration, file = paste0(i,"_Oncogenic_Signaling_Pathways_Alteration_dataset.rds"))
}                                          















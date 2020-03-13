#### 1.Generation_TCGA_Oncogenic_Signaling_Pathways_Alteration.R
# Paper:Oncogenic Signaling Pathways in The Cancer
### 1.Read Genomic Alteration Matrices_pathway_level table
Alteration_pathway_level <- read.csv("IM_Table S4. Genomic Alteration Matrices_pathway_level.csv")
Alteration_Alteration_level <- read.csv("Table_S4_Genomic_Alteration_Matrices_Figure3_Alteration_level.csv")
Alteration_gene_level <- read.csv("Table_S4_Genomic_Alteration_Matrices_Figure3_gene_level.csv")
### add uniform rownames
cbind(as.character(Alteration_pathway_level$SAMPLE_BARCODE),as.character(Alteration_Alteration_level$SAMPLE_BARCODE),
      as.character(Alteration_gene_level$SAMPLE_BARCODE))
sampleID <- as.character(Alteration_pathway_level$SAMPLE_BARCODE)
sampleID <- gsub("-",".",sampleID)
Alteration_pathway_level$rownames <- sampleID
Alteration_Alteration_level$rownames <- sampleID
Alteration_gene_level$rownames <- sampleID

rownames(Alteration_pathway_level)  <- sampleID
rownames(Alteration_Alteration_level)  <- sampleID
rownames(Alteration_gene_level)  <- sampleID
saveRDS(Alteration_pathway_level, "All__Oncogenic_Signaling_Alteration_pathway_level.rds")
#### 2.Read phenotype data #### 
phenoData <- read.csv("Table_S1_List_cancer type_and_subtype_Figure 1.csv") 
table(phenoData$DISEASE)
sampleID <- as.character(phenoData$SAMPLE_BARCODE)
sampleID <- gsub("-",".",sampleID)
rownames(phenoData) <- sampleID
phenoData$rownames <- sampleID
## How many cancer types
TCGA.cancer.types <- names(table(phenoData$DISEASE))
#### 3. Metadata generation ####
Alteration_pathway_level.metadata <- list(Paper = "Oncogenic Signaling Pathways in The Cancer",
                                          Table = "Table S4. Genomic Alteration Matrices_pathway_level.csv",
                                          doi =  "10.1016/j.cell.2018.03.035.")
#### 4.Generate_Oncogenic_Signaling_Pathways_Alteration_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- phenoData$DISEASE == i
  TCGA.pheno <- phenoData[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$rownames)
  ## Step3 Extract expression data by pheno data
  TCGA_df_pathway_level <- Alteration_pathway_level[TCGA.sampleID,]
  TCGA_df_Alteration_level <- Alteration_Alteration_level[TCGA.sampleID,]
  TCGA_df_gene_level <- Alteration_gene_level[TCGA.sampleID,]
  ## Step4 Biuld TCGA data sets
  TCGA_Oncogenic_Signaling_Alteration <- list(Alteration_pathway_level = TCGA_df,
                                              Alteration_Alteration_level = TCGA_df_Alteration_level,
                                              Alteration_gene_level = TCGA_df_gene_level,
                                              Pheno = TCGA.pheno,
                                              Alteration_pathway_level.metadata = Alteration_pathway_level.metadata)
  names(TCGA_Oncogenic_Signaling_Alteration)<-c(paste0(i,".Alteration.pathway.level"),
                                                paste0(i,".Alteration.Alteration.level"),
                                                paste0(i,".Alteration.gene.level"),
                                                "PhenoType",
                                                "Alteration_pathway_level.metadata")
  saveRDS(TCGA_Oncogenic_Signaling_Alteration, file = paste0(i,"_Oncogenic_Signaling_Alteration_dataset.rds"))
}                                          















#### 1.Generation_TCGA_Cancer_Aneuploidy_dataset.R
# Paper:Genomic and Functional Approaches to Understanding Cancer Aneuploidy
### 1.Read table
TCGA_Aneuploidy <- read.csv("Table_S2_Chromosome_Arm_Calls_and_Aneuploidy_Scores_Figure1.csv",
                                   header = TRUE)
sampleID <- as.character(TCGA_Aneuploidy$Sample)
sampleID <- gsub("-",".",sampleID)
TCGA_Aneuploidy$rownames <- sampleID
rownames(TCGA_Aneuploidy)  <- sampleID
### 2.Cancer Types
table(TCGA_Aneuploidy$Type)
TCGA.cancer.types <- names(table(TCGA_Aneuploidy$Type))
#### 3. Metadata generation ####
Aneuploidy.metadata <- list(Paper = "Genomic and Functional Approaches to Understanding Cancer Aneuploidy",
                            Table = "Table S2. Sample Chromosome Arm Calls and Aneuploidy Scores, Related to Figure 1",
                            doi =  "https://doi.org/10.1016/j.ccell.2018.03.007")
#### 4.Generate_Aneuploidy_dataset ####
#i="ACC"
for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- TCGA_Aneuploidy$Type == i
  TCGA_Aneuploidy_sub <- TCGA_Aneuploidy[TCGA_Index,]
  ## Step4 Biuld TCGA data sets
  TCGA_Aneuploidy_sub_list <- list(TCGA_Aneuploidy_sub = TCGA_Aneuploidy_sub,
                                              Aneuploidy.metadata = Aneuploidy.metadata)
  names(TCGA_Aneuploidy_sub_list)<-c(paste0(i,".Aneuploidy.score"),
                                                paste0(i,".Aneuploidy.metadata"))
  saveRDS(TCGA_Aneuploidy_sub_list, file = paste0(i,"_Aneuploidy_score_dataset.rds"))
}                                          


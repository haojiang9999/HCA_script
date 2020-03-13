#### 1_1_Generate_GI_Adenocarcinomas_Characteristics_dataset.R
### 1.Read table
## Paper: Comparative Molecular Analysis of Gastrointestinal Adenocarcinomas
GI_Adenocarcinomas_Characteristics <- read.csv("Table_S1_Summary_Table_of_Tumor_Sample_Characteristics.csv")
GI_Adenocarcinomas_Characteristics.metadata <- read.csv("Table_S1_Summary_Table_of_Tumor_Sample_Characteristics_metadata.csv")
# Add paper Info
paper <- data.frame("paper","Comparative Molecular Analysis of Gastrointestinal Adenocarcinomas")
colnames(paper) <- colnames(GI_Adenocarcinomas_Characteristics.metadata)
GI_Adenocarcinomas_Characteristics.metadata <- rbind(paper,GI_Adenocarcinomas_Characteristics.metadata)
# Read Epigenetic_Silencing_Calls
Epigenetic_Silencing_Calls <- read.csv("Table_S7_Epigenetic_Silencing_Calls_Figure_3.csv",row.names = 1)
Epigenetic_Silencing_Calls[1:5,1:5]
colnames(Epigenetic_Silencing_Calls)<- paste0("TCGA.",colnames(Epigenetic_Silencing_Calls))
Epigenetic_Silencing_Calls[1:5,1:5]
rownames(Epigenetic_Silencing_Calls)
## 2.Seperate the data
## How many cancer types
table(GI_Adenocarcinomas_Characteristics$TCGA.Project.Code)
TCGA.cancer.types <- names(table(GI_Adenocarcinomas_Characteristics$TCGA.Project.Code))

for(i in TCGA.cancer.types){
  # i = "COAD"
  ### Step1 separate data by cancer types ###
  TCGA_Index <- GI_Adenocarcinomas_Characteristics$TCGA.Project.Code == i
  GI_Adenocarcinomas_Characteristics_sub <- GI_Adenocarcinomas_Characteristics[TCGA_Index,]
  #### Step2 Convert and add rownames to the table #### 
  TCGA.sampleID <- as.character(GI_Adenocarcinomas_Characteristics_sub$SNP.Array.TCGA.Aliquot.ID)
  TCGA.sampleID <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1", TCGA.sampleID)
  TCGA.sampleID <- substr(TCGA.sampleID,1,nchar(TCGA.sampleID)-1)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  GI_Adenocarcinomas_Characteristics_sub$rownames <- TCGA.sampleID
  #### Step3 Extract Epigenetic_Silencing_Calls if possible
  #table(TCGA.sampleID %in% colnames(Epigenetic_Silencing_Calls))
  Epi.Samples <- TCGA.sampleID[TCGA.sampleID %in% colnames(Epigenetic_Silencing_Calls)]
  Epigenetic_Silencing_Calls_sub <- Epigenetic_Silencing_Calls[,Epi.Samples]
  ## Step4 Biuld TCGA data sets
  GI_Adenocarcinomas_Characteristics.list <- list(GI_Adenocarcinomas_Characteristics_sub = GI_Adenocarcinomas_Characteristics_sub,
                                                  GI_Adenocarcinomas_Characteristics.metadata = GI_Adenocarcinomas_Characteristics.metadata,
                                                  Epigenetic_Silencing_Calls_sub = Epigenetic_Silencing_Calls_sub)
  names(GI_Adenocarcinomas_Characteristics.list)<-c(paste0(i,".GI.Adenocarcinomas.Characteristics"),"GI.Adenocarcinomas.Characteristics.metadata",paste0(i,".Epigenetic.Silencing.Calls"))
  saveRDS(GI_Adenocarcinomas_Characteristics.list, file = paste0(i,"_GI_Adenocarcinomas_Characteristics_And_Epigenetic_Silencing_dataset.rds"))
}














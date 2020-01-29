##### Pheno Information construction
### Cell info addition
Pheno.merged<- Cor.Res.CV8000$Pheno.merged
Pheno.merged$CellType <- rep("No_data", nrow(Pheno.merged))
### Annotate the Normal and Tumor cells
# Normal ones
NMcellIndex <- grep("NC", rownames(Pheno.merged))
Pheno.merged$CellType[NMcellIndex] <- "Normal cells"
Pheno.merged[NMcellIndex, ]
NMcellIndex <- rownames(Pheno.merged) %in% colnames(NM.epithelial.COUNT.normal.uniGene)
Pheno.merged$CellType[NMcellIndex] <- "Normal cells"
Pheno.merged[NMcellIndex, ]
# Tumor ones
TcellIndex <- Pheno.merged$CellType == "No_data"
Pheno.merged$CellType[TcellIndex] <- "Tumor cells"
Pheno.merged <- Pheno.merged[,-1]
#### Complete annotation information
#names(Log2.expList)
#Tang.colon.cancer.FPKM.500.pheno
#Tang.colon.cancer.TPM.688.pheno
#NM.Epi.phenoType
#Tumor.Epi.phenoType
#### Tang info
Pheno.merged$patientID <- "Unknown"
Pheno.merged$sampleRegion <- "Unknown"
Pheno.merged[rownames(Tang.colon.cancer.FPKM.500.pheno),]$patientID <- as.character(Tang.colon.cancer.FPKM.500.pheno$patientID)
Pheno.merged[rownames(Tang.colon.cancer.TPM.688.pheno),]$patientID <- as.character(Tang.colon.cancer.TPM.688.pheno$patientID)
Pheno.merged[rownames(Tang.colon.cancer.FPKM.500.pheno),]$sampleRegion <- substr(as.character(Tang.colon.cancer.FPKM.500.pheno$sampleRegion)
                                                                         , start = 1, stop = 2)
Pheno.merged[rownames(Tang.colon.cancer.TPM.688.pheno),]$sampleRegion <- substr(as.character(Tang.colon.cancer.TPM.688.pheno$sampleRegion)
                                                                        , start = 1, stop = 2)
Pheno.merged$sampleRegion
#### 
Pheno.merged$Epi_cellTypes <- "Unknown"
NM.Epi.phenoType$Epi_cellTypes
Pheno.merged[rownames(NM.Epi.phenoType),]$Epi_cellTypes <- as.character(NM.Epi.phenoType$Epi_cellTypes)
Pheno.merged[rownames(Tumor.Epi.phenoType),]$Epi_cellTypes <- as.character(Tumor.Epi.phenoType$Epi_cellTypes)
View(Pheno.merged) # View the table

saveRDS(Pheno.merged, file = "scColon_Pheno.merged.rds")

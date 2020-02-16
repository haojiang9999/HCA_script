#### Read_all_COUNT_data.R
### 1.Read data
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE81861_human_colorectal_tumors/"
NM.all.COUNT <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_NM_all_cells_COUNT.csv.gz")))
Tumor.all.COUNT <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_tumor_all_cells_COUNT.csv.gz")))
NM.all.COUNT[1:5,1:5]
Tumor.all.COUNT[1:5,1:5]
cbind(NM.all.COUNT$X,Tumor.all.COUNT$X)
### 2.Gene annotation
GeneNames <- as.character(NM.all.COUNT$X)
GeneSymbol<-unlist(lapply(strsplit(GeneNames,"_"), '[[', 2))
Loc <- unlist(lapply(strsplit(GeneNames,"_"), '[[', 1))
# Some names not follow the pattern
GeneID <- unlist(lapply(strsplit(GeneNames,"_"), tail, 1))
GeneID2 <- gsub("\\..*","",GeneID)
GeneAnno <- data.frame(Loc= Loc, GeneSymbol = GeneSymbol,
                       GeneID = GeneID, GeneID2= GeneID2)
rownames(GeneAnno) <- GeneNames
### 3.Single cell phenoType
# Normal
NM.cells <- as.character(colnames(NM.all.COUNT[,-1]))
NM_ID <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 1))
All_cellTypes <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 2))
NM_color <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 3))
NM.all.phenoType <- data.frame(NM_ID = NM_ID, All_cellTypes = All_cellTypes,
                               NM_color = NM_color)
rownames(NM.all.phenoType) <- NM.cells
table(NM.all.phenoType$All_cellTypes)
# Tumor
tumor.cells <- as.character(colnames(Tumor.all.COUNT[,-1]))
tumor_ID <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 1))
All_cellTypes <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 2))
tumor_color <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 3))
tumor.all.phenoType <- data.frame(tumor_ID = tumor_ID, All_cellTypes = All_cellTypes,
                                  tumor_color = tumor_color)
rownames(tumor.all.phenoType) <- tumor.cells
table(tumor.all.phenoType$All_cellTypes)
### 4.Expression data process
# Normal
NM.all.COUNT.exp <- NM.all.COUNT[,-1] 
rownames(NM.all.COUNT.exp) <- GeneNames
NM.all.COUNT.exp[1:5,1:2]
# Tumor
Tumor.all.COUNT.exp <- Tumor.all.COUNT[,-1] 
rownames(Tumor.all.COUNT.exp) <- GeneNames
Tumor.all.COUNT.exp[1:5,1:2]
#### 5.Builed RCA epithelial COUNT dataset
RCA_all_COUNT_dataset <- list(NM.all.COUNT.exp = NM.all.COUNT.exp,
                              Tumor.all.COUNT.exp = Tumor.all.COUNT.exp,
                              NM.all.phenoType = NM.all.phenoType,
                              tumor.all.phenoType = tumor.all.phenoType,
                              GeneAnno = GeneAnno)
saveRDS(RCA_all_COUNT_dataset, file = "RCA_all_COUNT_dataset.rds")
#### 6.Normalize the COUNT data using the Seurat ####
library(Seurat)
# This methods like Tang TPM BUT for this data they do not uses UMI to 
# remove duplicate reads
NM.all.COUNT.normalized <- NormalizeData(NM.all.COUNT.exp,
                                          normalization.method = "RC",
                                            scale.factor = 1e6)
summary(Matrix::colSums(NM.all.COUNT.normalized))
Tumor.all.COUNT.normalized <- NormalizeData(Tumor.all.COUNT.exp,
                                               normalization.method = "RC",
                                               scale.factor = 1e6)
summary(Matrix::colSums(Tumor.all.COUNT.normalized))
RCA_all_COUNT_normalized_dataset <- list(NM.all.COUNT.normalized = NM.all.COUNT.normalized,
                                         Tumor.all.COUNT.normalized = Tumor.all.COUNT.normalized,
                                            NM.all.phenoType = NM.all.phenoType,
                                            tumor.all.phenoType = tumor.all.phenoType,
                                            GeneAnno = GeneAnno)
saveRDS(RCA_all_COUNT_normalized_dataset, file = "RCA_all_COUNT_normalized_dataset.rds")
####### Convert Ensembl ID to Gene symbol #######
## When multiple ID map to one symbol then select the highest expression ID as the symbol expression
rownames(NM.all.COUNT.normalized) <- GeneAnno$GeneID2
rownames(Tumor.all.COUNT.normalized) <- GeneAnno$GeneID2
NM.all.COUNT.normalized <- as.data.frame(as.matrix(NM.all.COUNT.normalized))
Tumor.all.COUNT.normalized <- as.data.frame(as.matrix(Tumor.all.COUNT.normalized))
# load the script
source("/data8t_4/JH/MyJobs/1_R_script/ExpMxID2Symbol.R")
NM.all.COUNT.normalized.uniGene <- ExpMxID2Symbol(NM.all.COUNT.normalized)
Tumor.all.COUNT.normalized.uniGene <- ExpMxID2Symbol(Tumor.all.COUNT.normalized)
#NM.epithelial.COUNT.normal.uniGene<- sapply(NM.epithelial.COUNT.normal.uniGene, unlist)
#NM.epithelial.COUNT.normal.uniGene<- as.data.frame(NM.epithelial.COUNT.normal.uniGene)
#Tumor.epithelial.COUNT.normal.uniGene<- sapply(Tumor.epithelial.COUNT.normal.uniGene, unlist)
#Tumor.epithelial.COUNT.normal.uniGene<- as.data.frame(Tumor.epithelial.COUNT.normal.uniGene)
# Build dataset 
RCA_all_COUNT_normalized_uniGene_dataset <- list(NM.all.COUNT.normalized.uniGene = NM.all.COUNT.normalized.uniGene,
                                                    Tumor.all.COUNT.normalized.uniGene = Tumor.all.COUNT.normalized.uniGene,
                                                    NM.all.phenoType = NM.all.phenoType,
                                                    tumor.all.phenoType = tumor.all.phenoType,
                                                    GeneAnno = GeneAnno)
saveRDS(RCA_all_COUNT_normalized_uniGene_dataset, file = "RCA_all_COUNT_normalized_uniGene_dataset.rds")


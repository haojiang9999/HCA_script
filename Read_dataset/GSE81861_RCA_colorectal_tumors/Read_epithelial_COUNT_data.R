### Read epithelial COUNT data
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE81861_human_colorectal_tumors/"
NM.epithelial.COUNT <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_NM_epithelial_cells_COUNT.csv.gz")))
tumor.epithelial.COUNT <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_tumor_epithelial_cells_COUNT.csv.gz")))
cbind(NM.epithelial.COUNT$X,tumor.epithelial.COUNT$X)
### Gene annotation
GeneNames <- as.character(NM.epithelial.COUNT$X)
GeneSymbol<-unlist(lapply(strsplit(GeneNames,"_"), '[[', 2))
Loc <- unlist(lapply(strsplit(GeneNames,"_"), '[[', 1))
# Some names not follow the pattern
GeneID <- unlist(lapply(strsplit(GeneNames,"_"), tail, 1))
GeneID2 <- gsub("\\..*","",GeneID)
GeneAnno <- data.frame(Loc= Loc, GeneSymbol = GeneSymbol,
                       GeneID = GeneID, GeneID2= GeneID2)
rownames(GeneAnno) <- GeneNames
### Single cell phenoType
# Normal
NM.cells <- as.character(colnames(NM.epithelial.COUNT[,-1]))
NM_ID <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 1))
Epi_cellTypes <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 2))
NM_color <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 3))
NM.Epi.phenoType <- data.frame(NM_ID = NM_ID, Epi_cellTypes = Epi_cellTypes,
                               NM_color = NM_color)
rownames(NM.Epi.phenoType) <- NM.cells
# Tumor
tumor.cells <- as.character(colnames(tumor.epithelial.COUNT[,-1]))
tumor_ID <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 1))
Epi_cellTypes <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 2))
tumor_color <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 3))
tumor.Epi.phenoType <- data.frame(tumor_ID = tumor_ID, Epi_cellTypes = Epi_cellTypes,
                                  tumor_color = tumor_color)
rownames(tumor.Epi.phenoType) <- tumor.cells
### Expression data process
# Normal
NM.epithelial.COUNT.exp <- NM.epithelial.COUNT[,-1] 
rownames(NM.epithelial.COUNT.exp) <- GeneNames
NM.epithelial.COUNT.exp[1:5,1:2]
# Tumor
tumor.epithelial.COUNT.exp <- tumor.epithelial.COUNT[,-1] 
rownames(tumor.epithelial.COUNT.exp) <- GeneNames
tumor.epithelial.COUNT.exp[1:5,1:2]
#### Builed RCA epithelial COUNT dataset
RCA_epithelial_COUNT_dataset <- list(NM.epithelial.COUNT.exp = NM.epithelial.COUNT.exp,
                                    tumor.epithelial.COUNT.exp = tumor.epithelial.COUNT.exp,
                                    NM.Epi.phenoType = NM.Epi.phenoType,
                                    tumor.Epi.phenoType = tumor.Epi.phenoType,
                                    GeneAnno = GeneAnno)
saveRDS(RCA_epithelial_COUNT_dataset, file = "RCA_epithelial_COUNT_dataset.rds")

#### Normalize the COUNT data using the Seurat ####
library(Seurat)
# Normalize expression data by Seurat using COUNT
library(Seurat)
# This methods like Tang TPM BUT for this data they do not uses UMI to 
# remove duplicate reads
NM.epithelial.COUNT.normal <- NormalizeData(NM.epithelial.COUNT.exp,
                                           normalization.method = "RC",
                                           scale.factor = 1e6)
summary(Matrix::colSums(NM.epithelial.COUNT.normal))

Tumor.epithelial.COUNT.normal <- NormalizeData(tumor.epithelial.COUNT.exp,
                                            normalization.method = "RC",
                                            scale.factor = 1e6)
summary(Matrix::colSums(Tumor.epithelial.COUNT.normal))

RCA_epithelial_COUNT_normal_dataset <- list(NM.epithelial.COUNT.normal = NM.epithelial.COUNT.normal,
                                     Tumor.epithelial.COUNT.normal = Tumor.epithelial.COUNT.normal,
                                     NM.Epi.phenoType = NM.Epi.phenoType,
                                     tumor.Epi.phenoType = tumor.Epi.phenoType,
                                     GeneAnno = GeneAnno)
saveRDS(RCA_epithelial_COUNT_normal_dataset, file = "RCA_epithelial_COUNT_normal_dataset.rds")

####### Convert Ensembl ID to Gene symbol #######
## When multiple ID map to one symbol then select the highest expression ID as the symbol expression
rownames(NM.epithelial.COUNT.normal) <- GeneAnno$GeneID2
rownames(Tumor.epithelial.COUNT.normal) <- GeneAnno$GeneID2
NM.epithelial.COUNT.normal <- as.data.frame(as.matrix(NM.epithelial.COUNT.normal))
Tumor.epithelial.COUNT.normal <- as.data.frame(as.matrix(Tumor.epithelial.COUNT.normal))
# load the script
source("/data8t_4/JH/MyJobs/1_R_script/ExpMxID2Symbol.R")
NM.epithelial.COUNT.normal.uniGene <- ExpMxID2Symbol(NM.epithelial.COUNT.normal)
Tumor.epithelial.COUNT.normal.uniGene <- ExpMxID2Symbol(Tumor.epithelial.COUNT.normal)
#NM.epithelial.COUNT.normal.uniGene<- sapply(NM.epithelial.COUNT.normal.uniGene, unlist)
#NM.epithelial.COUNT.normal.uniGene<- as.data.frame(NM.epithelial.COUNT.normal.uniGene)
#Tumor.epithelial.COUNT.normal.uniGene<- sapply(Tumor.epithelial.COUNT.normal.uniGene, unlist)
#Tumor.epithelial.COUNT.normal.uniGene<- as.data.frame(Tumor.epithelial.COUNT.normal.uniGene)
# Build dataset 
RCA_epithelial_COUNT_normal_uniGene_dataset <- list(NM.epithelial.COUNT.normal.uniGene = NM.epithelial.COUNT.normal.uniGene,
                                                    Tumor.epithelial.COUNT.normal.uniGene = Tumor.epithelial.COUNT.normal.uniGene,
                                                    NM.Epi.phenoType = NM.Epi.phenoType,
                                                    tumor.Epi.phenoType = tumor.Epi.phenoType,
                                                    GeneAnno = GeneAnno)
saveRDS(RCA_epithelial_COUNT_normal_uniGene_dataset, file = "RCA_epithelial_COUNT_normal_uniGene_dataset.rds")













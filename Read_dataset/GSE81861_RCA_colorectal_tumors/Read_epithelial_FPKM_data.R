### Read epithelial FPKM data
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE81861_human_colorectal_tumors/"
NM.epithelial.fpkm <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_NM_epithelial_cells_FPKM.csv.gz")))
tumor.epithelial.fpkm <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_tumor_epithelial_cells_FPKM.csv.gz")))
cbind(NM.epithelial.fpkm$X,tumor.epithelial.fpkm$X)
### Gene annotation
GeneNames <- as.character(NM.epithelial.fpkm$X)
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
NM.cells <- as.character(colnames(NM.epithelial.fpkm[,-1]))
NM_ID <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 1))
Epi_cellTypes <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 2))
NM_color <- unlist(lapply(strsplit(NM.cells,"__"), '[[', 3))
NM.Epi.phenoType <- data.frame(NM_ID = NM_ID, Epi_cellTypes = Epi_cellTypes,
                               NM_color = NM_color)
rownames(NM.Epi.phenoType) <- NM.cells
# Tumor
tumor.cells <- as.character(colnames(tumor.epithelial.fpkm[,-1]))
tumor_ID <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 1))
Epi_cellTypes <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 2))
tumor_color <- unlist(lapply(strsplit(tumor.cells,"__"), '[[', 3))
tumor.Epi.phenoType <- data.frame(tumor_ID = tumor_ID, Epi_cellTypes = Epi_cellTypes,
                               tumor_color = tumor_color)
rownames(tumor.Epi.phenoType) <- tumor.cells
### Expression data process
# Normal
NM.epithelial.fpkm.exp <- NM.epithelial.fpkm[,-1] 
rownames(NM.epithelial.fpkm.exp) <- GeneNames
NM.epithelial.fpkm.exp[1:5,1:2]
# Tumor
tumor.epithelial.fpkm.exp <- tumor.epithelial.fpkm[,-1] 
rownames(tumor.epithelial.fpkm.exp) <- GeneNames
tumor.epithelial.fpkm.exp[1:5,1:2]
#### Builed RCA epithelial FPKM dataset
RCA_epithelial_FPKM_dataset <- list(NM.epithelial.fpkm.exp = NM.epithelial.fpkm.exp,
                                    tumor.epithelial.fpkm.exp = tumor.epithelial.fpkm.exp,
                                    NM.Epi.phenoType = NM.Epi.phenoType,
                                    tumor.Epi.phenoType = tumor.Epi.phenoType,
                                    GeneAnno = GeneAnno)
saveRDS(RCA_epithelial_FPKM_dataset, file = "RCA_epithelial_FPKM_dataset.rds")




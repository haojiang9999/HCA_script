### read all_cells_FPKM
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE81861_human_colorectal_tumors/"
NM.all.fpkm <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_NM_all_cells_FPKM.csv.gz")))
tumor.all.fpkm <- read.csv(gzfile(paste0(filePath, "GSE81861_CRC_tumor_all_cells_FPKM.csv.gz")))
cbind(NM.epithelial.fpkm$X,tumor.epithelial.fpkm$X)
### Expression data process
# Normal
NM.all.fpkm.exp <- NM.all.fpkm[,-1] 
rownames(NM.all.fpkm.exp) <- GeneNames
NM.all.fpkm.exp[1:5,1:2]
# Tumor
tumor.all.fpkm.exp <- tumor.all.fpkm[,-1] 
rownames(tumor.all.fpkm.exp) <- GeneNames
tumor.all.fpkm.exp[1:5,1:2]
### Single cell phenoType
# Normal
NM.all.cells <- as.character(colnames(NM.all.fpkm[,-1]))
NM_ID <- unlist(lapply(strsplit(NM.all.cells,"__"), '[[', 1))
all_cellTypes <- unlist(lapply(strsplit(NM.all.cells,"__"), '[[', 2))
NM_color <- unlist(lapply(strsplit(NM.all.cells,"__"), '[[', 3))
NM.all.phenoType <- data.frame(NM_ID = NM_ID, all_cellTypes = all_cellTypes,
                               NM_color = NM_color)
table(all_cellTypes)
rownames(NM.all.phenoType) <- NM.all.cells
# Tumor
tumor.all.cells <- as.character(colnames(tumor.all.fpkm[,-1]))
tumor_ID <- unlist(lapply(strsplit(tumor.all.cells,"__"), '[[', 1))
all_cellTypes <- unlist(lapply(strsplit(tumor.all.cells,"__"), '[[', 2))
tumor_color <- unlist(lapply(strsplit(tumor.all.cells,"__"), '[[', 3))
tumor.all.phenoType <- data.frame(tumor_ID = tumor_ID, all_cellTypes = all_cellTypes,
                                  tumor_color = tumor_color)
rownames(tumor.all.phenoType) <- tumor.all.cells

#### save the datas
RCA_all_FPKM_dataset <- list(NM.all.fpkm.exp = NM.all.fpkm.exp,
                             tumor.all.fpkm.exp = tumor.all.fpkm.exp,
                             NM.all.phenoType = NM.all.phenoType,
                             tumor.all.phenoType = tumor.all.phenoType,
                                    GeneAnno = GeneAnno)
saveRDS(RCA_all_FPKM_dataset, file = "RCA_all_FPKM_dataset.rds")

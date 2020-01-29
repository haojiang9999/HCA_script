## Read Zhang Zemin Tcell in colon cancer
## Data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989
### 1.Read expression data TPM
TCell.S11138.TPM.file <- "/data8t_4/JH/GEO/GSE108989_T_cell_in_Colon_Cancer/GSE108989_CRC.TCell.S11138.TPM.txt.gz"
TCell.S11138.TPM <- read.table(TCell.S11138.TPM.file)
TCell.S11138.TPM[1:5,1:5]
### Extract gene expression data
TCell.S11138.TPM.exp <- TCell.S11138.TPM
# Change colnames
sampleNames <- as.matrix(TCell.S11138.TPM[1,])
colnames(TCell.S11138.TPM.exp) <- sampleNames
TCell.S11138.TPM.exp <- TCell.S11138.TPM.exp[-1,]# remove first row
TCell.S11138.TPM.exp[1:5,1:5] # Check
### Remove rows have no gene symbol 
## Which contain NA in symbol rows rwmove 23459 -> 23370 89 rows
TCell.S11138.TPM.exp <- TCell.S11138.TPM.exp[complete.cases(TCell.S11138.TPM.exp),]
## 2.Extract the gene annotation info
Gene.Anno <- TCell.S11138.TPM.exp[,1:2]
## Add rownames of gene symbol
# Finish the expression building
rownames(TCell.S11138.TPM.exp) <- TCell.S11138.TPM.exp$symbol
TCell.S11138.TPM.exp <- TCell.S11138.TPM.exp[,-(1:2)]
# convert factor to numeric
GSE108989_CRC_TCell.S11138.TPM.exp <- sapply(TCell.S11138.TPM.exp,function(x){
  as.numeric(levels(x))[x]
})
# add rownames
rownames(GSE108989_CRC_TCell.S11138.TPM.exp) <- Gene.Anno$symbol
GSE108989_CRC_TCell.S11138.TPM.exp[1:5,1:5]
GSE108989_CRC_TCell.S11138.TPM.exp <- as.data.frame(GSE108989_CRC_TCell.S11138.TPM.exp)
class(GSE108989_CRC_TCell.S11138.TPM.exp$`NP710-20180123`)
## 3.Original paper cluster and celltype annotation
cell_annotation <- read.table("/data8t_4/JH/GEO/GSE108989_T_cell_in_Colon_Cancer/GSE108989_cell_annotation.txt")
colnames(cell_annotation) <- as.matrix(cell_annotation[1,])
cell_annotation <- cell_annotation[-1,]
cell_annotation[1:5,]
rownames(cell_annotation) <- cell_annotation$UniqueCell_ID
## Read Cell Type annotation file
CellType_annotation <- read.csv("/data8t_4/JH/GEO/GSE108989_T_cell_in_Colon_Cancer/CellType_annotation.csv")
CellType_annotation <- CellType_annotation[,-c(3,4)]
cell_annotation <- dplyr::left_join(cell_annotation, CellType_annotation, by = "sampleType")
## Read cluster Info 
CellCluster_annotation <- read.csv("/data8t_4/JH/GEO/GSE108989_T_cell_in_Colon_Cancer/Zemin_GSE108989_Tcell_cluster_anno.csv")
cell_annotation <- dplyr::left_join(cell_annotation, CellCluster_annotation, by = "majorCluster")

### 4.Build GSE108989_CRC.TCell dataset
GSE108989_CRC.TCell.S11138.TPM.dataset <- list(GSE108989_CRC_TCell.S11138.TPM.exp = GSE108989_CRC_TCell.S11138.TPM.exp,
                                                Gene.Anno = Gene.Anno,
                                                Cell_annotation = cell_annotation)
saveRDS(GSE108989_CRC.TCell.S11138.TPM.dataset, file = "GSE108989_CRC.TCell.S11138.TPM.dataset.rds")














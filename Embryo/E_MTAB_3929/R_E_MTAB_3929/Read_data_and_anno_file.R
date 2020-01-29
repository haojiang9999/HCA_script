# loading data
counts <- read.table("../counts.txt")
FPKM <- read.table("../rpkm.txt")
# loading cell label from chen wenchang
sdrf <- read.delim2("../E-MTAB-3929.sdrf.txt")
head(sdrf)
load("../em_label1.RData")
load("../em_label2.RData")
labels <- label1$labels
lng <- label2$label.lng
head(lng[,1:10])
label2$label.lng[1,]
# create annotation file 
cellTypeNum = label2$label.lng[1,]
class(cellType)
cell.lineage <- cbind(cellName = colnames(FPKM), 
                      cellTypeNum = label2$label.lng[1,])
cell.lineage <- as.data.frame(cell.lineage)
cellType <- cbind(cellTypeNum = c("1","2","3","4"), 
      cellType = c("EPI","Prelineage","PE","TE"))
cellType <- as.data.frame(cellType)
cell.lineage <- dplyr::left_join(cell.lineage, cellType, by = "cellTypeNum")

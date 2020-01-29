##### Distance calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
Cor.Res.CV8000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.8000)
Cor.Res.CV4000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.4000)
Cor.Res.CV2000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.2000)
Cor.Res.CV1500 <- refCorMerge(Log2.expList, log2.scReference.list.CV.1500)
Cor.Res.CV1000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.1000)
Log2.Cor.Res.list <- list(Cor.Res.CV8000 = Cor.Res.CV8000,
     Cor.Res.CV4000 = Cor.Res.CV4000,
     Cor.Res.CV2000 = Cor.Res.CV2000,
     Cor.Res.CV1500 = Cor.Res.CV1500,
     Cor.Res.CV1000 = Cor.Res.CV1000,
     Pheno.merged = Pheno.merged)
saveRDS(Log2.Cor.Res.list, file = "Log2.Cor.Res.list.dataset.rds")

log2.CV1500.dataset <- list(Log2.expList = Log2.expList, 
                            log2.scReference.list.CV.1500 = log2.scReference.list.CV.1500)
saveRDS(log2.CV1500.dataset, file = "log2.CV1500.dataset.rds")
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



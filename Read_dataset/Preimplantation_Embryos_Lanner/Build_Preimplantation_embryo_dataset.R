#### Build Preimplantation embryo dataset
### 1.Read expression data download from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3929/
# read expression table
# loading data
Preim.embryo.COUNT <- read.table("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/counts.txt")
Preim.embryo.COUNT[1:10,1:10]
Preim.embryo.RPKM <- read.table("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/rpkm.txt")
summary(colSums(Preim.embryo.RPKM))
# Normalize expression data by Seurat using COUNT
library(Seurat)
# This methods like Tang TPM BUT for this data they do not uses UMI to 
# remove duplicate reads
Preim.embryo.COUNT.normal <- NormalizeData(Preim.embryo.COUNT,
                      normalization.method = "RC",
                      scale.factor = 1e6)
summary(Matrix::colSums(Preim.embryo.COUNT.normal))
### 2.Celltype annotation
sdrf <- read.delim2("/data8t_4/JH/MyJobs/Embryo/E_MTAB_3929/E-MTAB-3929.sdrf.txt")
table(sdrf$Characteristics.inferred.lineage.)
colnames(Preim.embryo.pheno)
rownames(sdrf) <- sdrf$Source.Name
Preim.embryo.pheno <- sdrf[colnames(Preim.embryo.COUNT.normal),]
Preim.embryo.pheno <- Preim.embryo.pheno[, c("Source.Name", "Comment.ENA_SAMPLE.",
                                             "Characteristics.individual." ,
                                             "Characteristics.developmental.stage.",
                                             "Characteristics.inferred.lineage.",
                                             "Characteristics.inferred.trophectoderm.subpopulation.",
                                             "Characteristics.inferred.pseudo.time." )]
# Cell type annotation
CellTypes <- gsub("not applicable", "Prelineage",Preim.embryo.pheno$Characteristics.inferred.lineage.)
Preim.embryo.pheno$CellTypes <- CellTypes
### 3.Gene annotation
Preim.embryo.geneNames <- rownames(Preim.embryo.COUNT.normal)
source("/data8t_4/JH/MyJobs/1_R_script/GeneSymbol2GeneID.R")
Preim.embryo.Gene.Anno <- GeneSymbol2GeneID(as.character(Preim.embryo.geneNames), toType = c("ENSEMBL","ENTREZID"))
length(unique(Preim.embryo.Gene.Anno$SYMBOL))
## 4.Build Adult colon dataset E-MTAB-3929
Pre_implantation_E_MTAB_3929_dataset <- list(Preim.embryo.exp.COUNT.normal = Preim.embryo.COUNT.normal, Preim.embryo.pheno = Preim.embryo.pheno,
                                             Preim.embryo.Gene.Anno = Preim.embryo.Gene.Anno)
#### Save the data set ####
saveRDS(Pre_implantation_E_MTAB_3929_dataset, file = "Pre_implantation_E_MTAB_3929_dataset.rds")

#### 1_Generation_MSigDB_Geneset_Object_V7_GMT.R
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.0_Release_Notes

library(GSEABase)
filePath <- "/stor/jianghao/Database/MSigDB/msigdb_v7.0_GMTs/"
### Read MSigDB metadata
MSigDB_V7_metadata <- read.csv("MSigDB_metadata.csv",header = F)
### 1.symbols.gmt
symbols.gmt <- list.files(filePath, pattern = "symbols.gmt")
files.symbols.gmt <- paste0(filePath,symbols.gmt)
symbols.gmt.list  <- sapply(files.symbols.gmt, function(x){
  getGmt(x)
})
names(symbols.gmt.list)  <- symbols.gmt
MSigDB_Geneset_Object_V7_GMT_symbols_dataset <- list(symbols.gmt.list = symbols.gmt.list, MSigDB_V7_metadata = MSigDB_V7_metadata)
saveRDS(MSigDB_Geneset_Object_V7_GMT_symbols_dataset, "MSigDB_gmt_V7_symbols_dataset.rds")

### 2.entryz.gmt
entrez.gmt <- list.files(filePath, pattern = "entrez.gmt")
files.entrez.gmt <- paste0(filePath,entrez.gmt)
files.entrez.gmt.list  <- sapply(files.entrez.gmt, function(x){
  getGmt(x)
})
names(files.entrez.gmt.list)  <- entrez.gmt
MSigDB_V7_entrez_dataset <- list(files.entrez.gmt.list = files.entrez.gmt.list, MSigDB_V7_metadata = MSigDB_V7_metadata)
saveRDS(MSigDB_V7_entrez_dataset, "MSigDB_gmt_V7_entrez_dataset.rds")




























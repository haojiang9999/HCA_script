#### 1_Generation_MSigDB_Geneset_Object_V7_GMT.R
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.0_Release_Notes

library(GSEABase)
filePath <- "/stor/jianghao/Database/MSigDB/msigdb_v7.0_GMTs/"
### 1.symbols.gmt
symbols.gmt <- list.files(filePath, pattern = "symbols.gmt")
files.symbols.gmt <- paste0(filePath,symbols.gmt)
symbols.gmt.list  <- sapply(files.symbols.gmt, function(x){
  getGmt(x)
})
names(symbols.gmt.list)  <- symbols.gmt
saveRDS(symbols.gmt.list, "MSigDB_Geneset_Object_V7_GMT_symbols.rds")

### 2.
entrez.gmt <- list.files(filePath, pattern = "entrez.gmt")
files.entrez.gmt <- paste0(filePath,entrez.gmt)
files.entrez.gmt.list  <- sapply(files.entrez.gmt, function(x){
  getGmt(x)
})
names(files.entrez.gmt.list)  <- entrez.gmt
saveRDS(files.entrez.gmt.list, "MSigDB_Geneset_Object_V7_GMT_entrez.rds")




























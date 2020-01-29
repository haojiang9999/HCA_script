#lOAD DATA FOR FUNs
#### 1.WikiPathways analysis ####
wpgmtfile <- "Database/WikiPathways/wikipathways-20190810-gmt-Homo_sapiens.gmt"
#### 2.Cell Marker analysis ####
cellMarkerFile <- "/stor/jianghao/Database/Cell_Marker/Human_cell_markers.txt"
#### 3.MSigDb analysis ####
#The MSigDB is a collection of annotated gene sets, 
#it include 8 major collections:
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures
## database path

msigdbFilePath <- "Database/MSigDB/"
# load MSigDB
entrez.gmt <- list.files(msigdbFilePath, pattern = ".entrez.gmt")
msigdbFiles <- file.path(msigdbFilePath, entrez.gmt)
## Loading database
#c1 <- read.gmt(msigdbFiles[1])
c2 <- read.gmt(msigdbFiles[2])
c3 <- read.gmt(msigdbFiles[3])
c4 <- read.gmt(msigdbFiles[4])
#c5 <- read.gmt(msigdbFiles[5])
c6 <- read.gmt(msigdbFiles[6])
c7 <- read.gmt(msigdbFiles[7])
h <- read.gmt(msigdbFiles[8])
# 

#### 4.Gene Ontology Analysis ####
#### JH version enrichGO prepare
GO_DATA_BP <- readRDS("Database/GO_DATA/GO_DATA_BP.rds")
GO_DATA_CC <- readRDS("Database/GO_DATA/GO_DATA_CC.rds")
GO_DATA_MF <- readRDS("Database/GO_DATA/GO_DATA_MF.rds")
### load clusterProfiler_JH self made code
library(devtools)
packPath <- "clusterProfiler_JH/"
load_all(packPath)







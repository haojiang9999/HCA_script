#### 1.WikiPathways analysis####
## test input data 
csvFile <- read.csv("test4.csv")

x <- unlist(csvFile[,2:5])
input <- as.character(unique(x))
# convert geneSymble to geneID
inputID <- geneSymbol2GeneID(input)
inputID <- inputID$ENTREZID
## start wiki annotation
library(magrittr)
library(clusterProfiler)
# read wiki gmt files
wpgmtfile <- "/stor/jianghao/Database/WikiPathways/wikipathways-20190810-gmt-Homo_sapiens.gmt"
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
ewp <- enricher(inputID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
# add background genes!!!!
head(ewp)

#### 2.Cell Marker analysis ####
cellMarkerFile <- "/stor/jianghao/Database/Cell_Marker/Human_cell_markers.txt"
cell_markers <- vroom::vroom(cellMarkerFile) %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers

y <- enricher(inputID, TERM2GENE=cell_markers, minGSSize=1)
DT::datatable(as.data.frame(y))

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
msigdbFilePath <- "/stor/jianghao/Database/MSigDB/"
entrez.gmt <- list.files(msigdbFilePath, pattern = ".entrez.gmt")
msigdbFiles <- file.path(msigdbFilePath, entrez.gmt)
c1 <- read.gmt(msigdbFiles[1])
c2 <- read.gmt(msigdbFiles[2])
c3 <- read.gmt(msigdbFiles[3])
c4 <- read.gmt(msigdbFiles[4])
c5 <- read.gmt(msigdbFiles[5])
h <- read.gmt(msigdbFiles[6])

egmtc1 <- enricher(gene, TERM2GENE=c1)
egmtc2 <- enricher(gene, TERM2GENE=c2)
egmtc3 <- enricher(gene, TERM2GENE=c3)
egmtc4 <- enricher(gene, TERM2GENE=c4)
egmtc5 <- enricher(gene, TERM2GENE=c5)
egmth <- enricher(gene, TERM2GENE=h)

msigdb2gene <- read.gmt(msigdbFile)
em <- enricher(inputID, TERM2GENE=msigdb2gene)
head(em)
egmtc1@result
head(egmtc5@result)






#### 4.Gene Ontology Analysis ####
egoCC <- enrichGO(gene          = inputID,
                #universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(egoCC)

egoMF <- enrichGO(gene          = inputID,
                  #universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoMF)

egoBP <- enrichGO(gene          = inputID,
                  #universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoBP)
test <- dropGO(egoBP, level = 6)
test@result
head(test@result)


#### 5. KEGG analysis ####
library(clusterProfiler)
search_kegg_organism('hsa', by='kegg_code') # hsa for human
HomoSapiens <- search_kegg_organism('Homo sapiens', by='scientific_name')
dim(HomoSapiens)
head(HomoSapiens)

# KEGG over-representation test
inputID
gene <- names(geneList)[abs(geneList) > 2]
# kegg_pATHWAY
kk <- enrichKEGG(gene         = inputID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05 )
head(kk)
# KEGG Module over-representation test
mkk <- enrichMKEGG(gene = inputID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
head(mkk)
mkk@result

#### 6.Reactome pathway analysis####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")
library(ReactomePA)

ere<- enrichPathway(inputID)
head(ere)
enrichPathway


ere@result

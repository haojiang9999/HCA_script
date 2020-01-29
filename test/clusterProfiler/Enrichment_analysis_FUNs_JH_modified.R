##### Enrichment analysis using clusterProfiler #######
# test input data
# convert geneSymble to geneID
inputID <- geneSymbol2GeneID(input)
inputID <- inputID$ENTREZID

#### 1.WikiPathways analysis ####
wpgmtfile <- "/stor/jianghao/Database/WikiPathways/wikipathways-20190810-gmt-Homo_sapiens.gmt"

WikiEnrich <- function(inputID, wpgmtfile= wpgmtfile, background ){
  wp2gene <- read.gmt(wpgmtfile)
  wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  ewp <- enricher(inputID, TERM2GENE = wpid2gene, 
                  TERM2NAME = wpid2name,
                  universe = background)
  return(ewp@result)
}
#### 2.Cell Marker analysis ####
cellMarkerFile <- "/stor/jianghao/Database/Cell_Marker/Human_cell_markers.txt"
#
CellMarkerEnrich <- function(inputID ,cellMarkerFile = cellMarkerFile,
                             background # background genes
){
  if (!require("DT")) install.packages('vroom')
  cell_markers <- vroom::vroom(cellMarkerFile) %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
  #cell_markers
  ecm <- enricher(inputID, TERM2GENE=cell_markers, minGSSize=1,
                  universe = background )
  return(ecm@result)
  
}
# DT::datatable(as.data.frame( ))

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

MSigDBEnrich <- function(inputID,msigdbFilePath = msigdbFilePath ,background){
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
  MSigDBresaults <- list()
  #egmtc1 <- enricher(gene = inputID, universe = background,TERM2GENE=c1)
  egmtc2 <- enricher(gene = inputID, universe = background,TERM2GENE=c2)
  egmtc3 <- enricher(gene = inputID, universe = background, TERM2GENE=c3)
  egmtc4 <- enricher(gene = inputID, universe = background, TERM2GENE=c4)
  #egmtc5 <- enricher(gene = inputID, universe = background, TERM2GENE=c5)
  egmtc6 <- enricher(gene = inputID, universe = background, TERM2GENE=c6)
  egmtc7 <- enricher(gene = inputID, universe = background, TERM2GENE=c7)
  egmth <- enricher(gene = inputID, universe = background, TERM2GENE=h)
  #MSigDBresaults$C1_positional_gene_sets <- egmtc1@result
  MSigDBresaults$C2_curated_gene_sets <- egmtc2@result
  MSigDBresaults$C3_motif_gene_sets <- egmtc3@result
  MSigDBresaults$C4_computational_gene_sets <- egmtc4@result
  #MSigDBresaults$C5_GO_gene_sets <- egmtc5@result
  MSigDBresaults$C6_oncogenic_signatures <- egmtc6@result
  MSigDBresaults$C7_immunologic_signatures <- egmtc7@result
  MSigDBresaults$H_hallmark_gene_sets <- egmth@result
  return(MSigDBresaults) 
}


#### 4.Gene Ontology Analysis ####
#### JH version enrichGO prepare
GO_DATA_BP <- readRDS("/stor/jianghao/Database/GO_DATA/GO_DATA_BP.rds")
GO_DATA_CC <- readRDS("/stor/jianghao/Database/GO_DATA/GO_DATA_CC.rds")
GO_DATA_MF <- readRDS("/stor/jianghao/Database/GO_DATA/GO_DATA_MF.rds")
### load clusterProfiler_JH self made code
library(devtools)
packPath <- "/data8t_4/JH/MyJobs/test/R_package_process/clusterProfiler_JH"
load_all(packPath)

GOEnrich <-  function(inputID, background, ...){
  
  GOresault <- list()
  ## cellular component (CC)
  egoCC <- enrichGO_JH(gene = inputID,
                    universe      = background,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  ## molecular function (MF)
  egoMF <- enrichGO_JH(gene = inputID,
                    universe      = background,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  ## biological process (BP)
  egoBP <- enrichGO_JH(gene = inputID,
                    universe      = background,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  ## drop GO terms duplicated or at specific level
  egoCC <- dropGO(egoCC)
  egoMF <- dropGO(egoMF)
  egoBP <- dropGO(egoBP)
  GOresault$egoCC <- egoCC@result
  GOresault$egoMF <- egoMF@result
  GOresault$egoBP <- egoBP@result
  return(GOresault)
}
#### 5.KEGG analysis ####
# kegg_pATHWAY

KEGGEnrich <- function(inputID, background){
  KEGGresaults <- list()
  # KEGG Pathway
  kk <- enrichKEGG(gene         = inputID,
                   universe = background,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  # KEGG Module over-representation test
  mkk <- enrichMKEGG(gene = inputID,
                     universe = background,
                     organism = 'hsa',
                     pvalueCutoff = 0.05)
  # output resault
  KEGGresaults$KEGG_Pathway <- kk@result
  KEGGresaults$KEGG_MODULE <- mkk@result
  return(KEGGresaults)
}

#### 6.Reactome pathway analysis ####
ReactomeEnrich <- function(inputID, background){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("ReactomePA"))
    BiocManager::install("ReactomePA")
  ere <- ReactomePA::enrichPathway(gene = inputID,
                                   universe = background,
                                   pvalueCutoff = 0.05)
  return(ere@result)
}








#### convert gene symbol to geneID or others ####
#input = geneNamesFPKM
GeneSymbol2GeneID <- function(input, fromType="SYMBOL", 
                              toType="ENTREZID", 
                              OrgDb="org.Hs.eg.db"){
  ### depending packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("org.Hs.eg.db")) 
    BiocManager::install("org.Hs.eg.db")
  if (!require("clusterProfiler")) 
    BiocManager::install("clusterProfiler")
  require(clusterProfiler)
  require(org.Hs.eg.db)
  #### Step1. convert geneSymbol to geneID
  ids1 = bitr(input, fromType=fromType, toType=toType, 
              OrgDb=OrgDb,drop = F)
 
  #### Step2. Find geneSymbols not recognized in step1
  leftSymbol <- ids1[is.na(ids1[,2]),][,1]
 
  # check whethere all geneSymbol had been finded
  
  if (length(leftSymbol) > 0) {
    ids2 = bitr(leftSymbol, fromType="ALIAS", toType=toType, OrgDb=OrgDb,drop = FALSE)
    colnames(ids2)[1] <- "SYMBOL"
    ids <- rbind(ids1,ids2)
    
  } else {ids <- ids1}
  ids <- na.omit(ids) # remove na row
  return(ids)
}

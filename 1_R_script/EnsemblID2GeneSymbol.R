#### convert ENSEMBL ID to gene symbol ####
#input = geneNamesFPKM
EnsemblID2GeneSymbol <- function(input, fromType="ensembl_gene_id", 
                              toType="hgnc_symbol"){
  ### depending packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("biomaRt")) 
    BiocManager::install("biomaRt")
  require(biomaRt)
  #### Step1.Loading dataset from biomaRt
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  e2s <- getBM(attributes=c(fromType, toType),
               filters= fromType,
               values= input,
               mart=ensembl)
  return(e2s)
  
}


#### Function convert expression matrix rowname Ensembl ID to Gene symbol
ExpMxID2Symbol <- function(exp.tb){
  ENSEMBL <- rownames(exp.tb)
  source("/data8t_4/JH/MyJobs/1_R_script/EnsemblID2GeneSymbol.R")
  eg <- EnsemblID2GeneSymbol(ENSEMBL)
  ## remove EnsemblID not have Symbol
  eg.sub <- eg[nchar(eg$hgnc_symbol) != 0,]
  # remove same ensembl map to different symbol
  eg.sub <- eg.sub[!duplicated(eg.sub$ensembl_gene_id),]
  rownames(eg.sub) <- eg.sub$ensembl_gene_id
  ## 1.select ENSEMBL had a SYMBOL paire
  ENSEMBL.dup <- eg.sub[duplicated(eg.sub$hgnc_symbol),]$ensembl_gene_id
  ### Step1.Build Non-dup symbol expression matrix
  # Symbol of dup or non-dup genes
  Dup.genes <- as.character(unique(eg.sub[duplicated(eg.sub$hgnc_symbol),]$hgnc_symbol))
  #NonDup.genes <-as.character(unique(GeneAnno.sub[!duplicated(GeneAnno.sub$GeneSymbol),]$GeneSymbol))
  NonDup.genes <-as.character(eg.sub[!eg.sub$hgnc_symbol %in% Dup.genes,]$hgnc_symbol)
  ### Expression dataframe of Non-dup genes
  Index.NonDup<- eg.sub[!eg.sub$hgnc_symbol %in% Dup.genes,]$ensembl_gene_id## remove dup genes
  table(duplicated(Index.NonDup))
  exp.tb.uni1 <- exp.tb[Index.NonDup,]
  ## change ensembl_gene_id to hgnc_symbol
  rownames(exp.tb.uni1) <- eg.sub[Index.NonDup,]$hgnc_symbol
  ### Step2.Build duplicated ensembl_gene_id expression matrix using Max expression one
  exp.tb.uni2 <- sapply(Dup.genes, FUN = function(symbol){
    #symbol= "RNU6-545P"
    Index.Dup<- rownames(eg.sub[eg.sub$hgnc_symbol %in% symbol,])
    exp.uni2 <- exp.tb[Index.Dup,]
    # Find Ensembl ID that have highest sum value
    exp.rowsum <- rowSums(exp.uni2)
    Name.Max <- names(exp.rowsum[exp.rowsum == max(exp.rowsum)][1]) # when all was 0 select only 1 ID
    exp.uni2 <- exp.uni2[Name.Max,]
    exp.uni2 <- as.data.frame(exp.uni2)
  })
  exp.tb.uni2 <- t(exp.tb.uni2)
  #### Step3.Combine two unique table
  exp.tb.uniGene <- rbind(exp.tb.uni1,exp.tb.uni2)
  geneNames <- rownames(exp.tb.uniGene)
  exp.tb.uniGene <- sapply(exp.tb.uniGene, unlist)
  exp.tb.uniGene<- as.data.frame(exp.tb.uniGene)
  rownames(exp.tb.uniGene) <- geneNames
  return(exp.tb.uniGene)
  
}

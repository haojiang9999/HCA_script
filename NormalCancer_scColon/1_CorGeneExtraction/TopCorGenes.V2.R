#### Step4 Find top N genes with highest CorValue in each cell types ##### 
#CorValueResaults <- CorValueByGenes_V3(inputCells, Filtered_referPanel_of_mix)
### construc function
TopCorGenes <- function(CorValueResaults, # input list from CorValueByGenes_V3
                        TopGeneNum = 200  # Number of top correlated genes in each cell types
){

  #### Find top TopGeneNum genes with highest CorValue in each cell types ####  
  TopCorGenes_mx  <- lapply(CorValueResaults, function(cellGeneCor){
    # Gene names of the df
    geneNames <- rownames(cellGeneCor)
    # select top genes in different cell types
    apply(cellGeneCor,2, function(x){
      TopIndex <- order(x, decreasing = T)[1:TopGeneNum] # how many top genes
      geneNames[TopIndex]
    })
  })

  return(TopCorGenes_mx)
}

#TopNcellTypes <- TopNcellTypes_df[[1]]


#### Generate a table of top genes and top cell types from sorted genes ####
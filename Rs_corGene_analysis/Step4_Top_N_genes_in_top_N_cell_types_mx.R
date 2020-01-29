#### Step4 Find top ** genes with highest CorValue in top * cell types ##### 
# Load test data
inputCells = GSE97693_Tang_TPM_cells[,1:20] # drop very important
referencePanel = referPanel_of_mix
CorValueResaults <- CorValueByGenes_V3(inputCells, Filtered_referPanel_of_mix)
### construc function
TopCorGenes <- function(CorValueResaults, # input list from CorValueByGenes_V3
                        cellTypeNum = 5,  # Number of top correlated cell types
                        TopGeneNum = 200  # Number of top correlated genes in each cell types
                        ){
  #### step1 select high correlation score reference cell types ####
  ## calculate correlation score of each cell through deferent reference cell types
  cellCor <- lapply(CorValueResaults, FUN = colSums)
  TopNcellTypes <- lapply(cellCor, function(x){
    names(head(sort(x, decreasing = T),cellTypeNum))
  })
  #### step2 select top N reference cell types related genes ####
  TopNcellTypes_df <- list()
  for(i in 1:length(CorValueResaults)){
    TopNcellTypes_df[[i]] <- CorValueResaults[[i]][,TopNcellTypes[[i]]]
  }
  names(TopNcellTypes_df) <- names(CorValueResaults)
  #### step3 find top ** genes with highest CorValue in top 5 cell types ####  
  TopCorGenes_mx  <- lapply(TopNcellTypes_df, function(TopNcellTypes){
    # Gene names of the df
    geneNames <- rownames(TopNcellTypes)
    # select top genes in different cell types
    apply(TopNcellTypes,2, function(x){
      TopIndex <- order(x, decreasing = T)[1:TopGeneNum] # how many top genes
      geneNames[TopIndex]
    })
  })
  return(TopCorGenes_mx)
}

#TopNcellTypes <- TopNcellTypes_df[[1]]


#### Generate a table of top genes and top cell types from sorted genes ####

                        
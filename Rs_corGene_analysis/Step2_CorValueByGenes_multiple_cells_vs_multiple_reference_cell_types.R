##### Step2 calculate CorValueByGenes multiple cells to multiple reference cell types ######
## input for test
inputCells = GSE97693_Tang_TPM_cells[,1:20] # drop very important
referencePanel = referPanel_of_mix[1:10]
i=1

#### fuction construction
######### Full version but output too large #########
CorValueByGenes <- function(inputCells, referencePanel){
  singleCellCorValue <- list()
  ### calculate CorValue by cells
  for(i in 1:length(colnames(inputCells))){
    inputCell = inputCells[,i,drop = F] # drop = F to keep rownames
    singleCellCorValue[[i]] <- FindCorValueByGenes(inputCell, referencePanel)
  }
  
}
######### only use single cell > 0 genes version V3 FindCorValueByGenes #########
CorValueByGenes_V3 <- function(inputCells, referencePanel){
  singleCellCorValue <- list()
  ### calculate CorValue by cells
  for(i in 1:length(colnames(inputCells))){
    inputCell = inputCells[,i,drop = F] # drop = F to keep rownames
    singleCellCorValue_output <- FindCorValueByGenes(inputCell, referencePanel)
    ### select CellCorValue from singleAbove0 version 
    singleCellCorValue[[i]] <- singleCellCorValue_output$singleAbove0_df
  }
  names(singleCellCorValue) <- colnames(inputCells)
  return(singleCellCorValue)
  
}










####### Step1 calculate CorValueByGenes 1 cell to multiple reference cell types #######
## input for test
#inputCell = GSE97693_Tang_TPM_cells[,1,drop = F] # drop = F to keep rownames
#referencePanel = referPanel_of_mix
#i=1

# fuction construction
FindCorValueByGenes <- function(inputCell, referencePanel){
  ### pre step: find interscet genes
  geneList <- intersect(rownames(inputCell),rownames(referencePanel))
  inputCell <- inputCell[geneList,]
  referencePanel <- referencePanel[geneList,]
  corValueGenes <- list() # list to contain resaults
  ############### Version: reference and single cell all genes #####################
  for(i in 1:length(colnames(referencePanel))) {
    covForGenes <-  ((inputCell-mean(inputCell))*(referencePanel[,i]-mean(referencePanel[,i]))/length(referencePanel[,i]))
    sdForGenePairs <- (sd(nonZeroInputCell)*sd(referencePanel[,i]))
    corValueGenes[[i]] <- covForGenes/sdForGenePairs
    # add reference cell type names
    names(corValueGenes[[i]]) <- rownames(referencePanel)
   
  }
  names(corValueGenes) <- colnames(referencePanel) # add reference cell type names
  corValueGenes_V1 <- corValueGenes                # how to calculate the corValueGenes
  ############### Version: reference + single cell > 0 genes #######################
  for(i in 1:length(colnames(referencePanel))) {
    nonZer0Index <- inputCell + referencePanel[,i] > 0 # find non-zero in both vectors
    geneList.2 <- geneList[nonZer0Index]
    nonZeroInputCell <- inputCell[nonZer0Index]
    nonZeroReference <- referencePanel[,i][nonZer0Index]
    covForGenes <-  ((nonZeroInputCell-mean(nonZeroInputCell))*(nonZeroReference-mean(nonZeroReference))/length(nonZeroReference))
    sdForGenePairs <- (sd(nonZeroInputCell)*sd(nonZeroReference))
    corValueGenes[[i]] <- covForGenes/sdForGenePairs
    # add gene names to corValueGenes
    names(corValueGenes[[i]]) <- geneList.2
  }
  names(corValueGenes) <- colnames(referencePanel)
  corValueGenes_V2 <- corValueGenes
  ############### Version: Gene expression of single cell > 0 genes #######################
  for(i in 1:length(colnames(referencePanel))) {
    nonZer0Index <- inputCell > 0                # single cell > 0 genes
    nonZeroInputCell <- inputCell[nonZer0Index]
    nonZeroReference <- referencePanel[,i][nonZer0Index]
    geneList.2 <- geneList[nonZer0Index]
    # calculate corValueGenes
    covForGenes <-  ((nonZeroInputCell-mean(nonZeroInputCell))*(nonZeroReference-mean(nonZeroReference))/length(nonZeroReference))
    sdForGenePairs <- (sd(nonZeroInputCell)*sd(nonZeroReference))
    corValueGenes[[i]] <- covForGenes/sdForGenePairs
    # add gene names to corValueGenes
    names(corValueGenes[[i]]) <- geneList.2
  }
  names(corValueGenes) <- colnames(referencePanel)
  corValueGenes_V3 <- corValueGenes
  ### make corValueGenes_V3 a matrix
  corValueGenes_V3_df <- do.call("cbind",corValueGenes_V3)
  
  ### combine 3 version resaults
  corValueGenesDataSet <- list(allGenes = corValueGenes_V1, BothAbove0 = corValueGenes_V2, 
                               singleAbove0 = corValueGenes_V3, singleAbove0_df = corValueGenes_V3_df)
  return(corValueGenesDataSet)
}












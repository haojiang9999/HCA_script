###### Normal Cancer core script #######
### Function calculate corelation Matrix
corMatrix <- function(input, reference, method = "pearson"){
  geneList <- intersect(rownames(input),rownames(reference))
  df_combine = cbind(input[geneList,],reference[geneList,])
  cor_matrix = cor(df_combine,method = method)
  cor_matrix = as.data.frame(cor_matrix)
  cell_type_matrix = cor_matrix[(dim(input)[2]+1):dim(cor_matrix)[2],1:dim(input)[2]]
}

####### Step1 calculate CorValueByGenes 1 cell to multiple reference cell types #######
FindCorValueByGenes <- function(inputCell, referencePanel){
  ### pre step: find interscet genes
  geneList <- intersect(rownames(inputCell),rownames(referencePanel))
  inputCell <- inputCell[geneList,]
  referencePanel <- referencePanel[geneList,]
  corValueGenes <- list() # list to contain resaults
  ############### Version: reference and single cell all genes #####################
  # i=1
  for(i in 1:length(colnames(referencePanel))) {
    covForGenes <-  ((inputCell-mean(inputCell))*(referencePanel[,i]-mean(referencePanel[,i]))/length(referencePanel[,i]))
    sdForGenePairs <- (sd(inputCell)*sd(referencePanel[,i]))
    corValueGenes[[i]] <- covForGenes/sdForGenePairs
    # add reference cell type names
    names(corValueGenes[[i]]) <- colnames(referencePanel)
    
  }
  names(corValueGenes) <- colnames(referencePanel) # add reference cell type names
  corValueGenes_V1 <- corValueGenes                # how to calculate the corValueGenes
  return(corValueGenes_V1)
}
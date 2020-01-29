### calculate CorValueByGenes 1 cell to multiple reference cell types

## Load sample
inputCell = GSE97693_Tang_TPM_cells[,1,drop = F]
referencePanel = Filtered_referPanel_of_mix
# find interscet genes
geneList <- intersect(rownames(inputCell),rownames(referencePanel))
i=1
inputCell <- inputCell[geneList,]
referencePanel <- referencePanel[geneList,]
corValueGenes <- list()

############### Version: reference and single cell all genes #####################
for(i in 1:length(colnames(referencePanel))) {
  covForGenes <-  ((inputCell-mean(inputCell))*(referencePanel[,i]-mean(referencePanel[,i]))/length(referencePanel[,i]))
  sdForGenePairs <- (sd(nonZeroInputCell)*sd(referencePanel[,i]))
  corValueGenes[[i]] <- covForGenes/sdForGenePairs
  #print(i)
}
############### Version: reference + single cell > 0 genes #######################
for(i in 1:length(colnames(referencePanel))) {
  nonZer0Index <- inputCell + referencePanel[,i] > 0 # find non-zero in both vectors
  nonZeroInputCell <- inputCell[nonZer0Index]
  nonZeroReference <- referencePanel[,i][nonZer0Index]
  covForGenes <-  ((nonZeroInputCell-mean(nonZeroInputCell))*(nonZeroReference-mean(nonZeroReference))/length(nonZeroReference))
  sdForGenePairs <- (sd(nonZeroInputCell)*sd(nonZeroReference))
  corValueGenes[[i]] <- covForGenes/sdForGenePairs
  #print(i)
}
names(corValueGenes) <- colnames(referencePanel)
test<- do.call("cbind", corValueGenes)
dim(test)
lapply(corValueGenes, length)
cor(inputCell,referencePanel)
summary(cor(inputCell,referencePanel))
plot(corValueGenes[[1]])
plot(sort(corValueGenes[[1]], decreasing = T))
summary(corValueGenes[[1]])
length(corValueGenes[[1]])

hist(corValueGenes[[1]])
corValueGenes[[1]][corValueGenes[[1]] <0 ]

############### Version: Gene expression of single cell > 0 genes #######################
for(i in 1:length(colnames(referencePanel))) {
  nonZer0Index <- inputCell + referencePanel[,i] > 0 # find non-zero in both vectors
  nonZeroInputCell <- inputCell[nonZer0Index]
  nonZeroReference <- referencePanel[,i][nonZer0Index]
  geneList.2 <- geneList[nonZer0Index]
  # Non-Zero in single cells
  nonZer0SingleIndex <- nonZeroInputCell > 0
  geneList.3 <- geneList.2[nonZer0SingleIndex]
  nonZer0SingleIndex_InputCell <- nonZeroInputCell[nonZer0SingleIndex]
  nonZer0SingleIndex_referencePanel <- nonZeroReference[nonZer0SingleIndex]
  covForGenes <-  ((nonZer0SingleIndex_InputCell-mean(nonZer0SingleIndex_InputCell))*(nonZer0SingleIndex_referencePanel-mean(nonZer0SingleIndex_referencePanel))/length(nonZer0SingleIndex_referencePanel))
  sdForGenePairs <- (sd(nonZer0SingleIndex_InputCell)*sd(nonZer0SingleIndex_referencePanel))
  corValueGenes[[i]] <- covForGenes/sdForGenePairs
  #print(i)
}
corValueGenes[[i]] 
summary(corValueGenes[[i]] )
plot(sort(corValueGenes[[1]], decreasing = T)[1:200])
abline(h = 0.0001)
summary(corValueGenes[[1]] > 0)
summary(sort(corValueGenes[[1]], decreasing = T)[1:1000])
test <- geneList.3[order(corValueGenes[[1]], decreasing = T)][1:100]
write.csv(test, file = "test.csv")
colnames(referencePanel)

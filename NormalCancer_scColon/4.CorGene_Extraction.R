#### Correlation gene extraction 
test1 <- CorValueByGenes(mxList[["Tang.colon.cancer.FPKM.500"]],ReferenceList[["Tang.Adult.colon"]])
inputCell = mxList[["Tang.colon.cancer.FPKM.500"]]
referencePanel = ReferenceList[["Tang.Adult.colon"]]
test1 <- CorValueByGenes(inputCell,referencePanel) 
test2 <- FindCorValueByGenes(inputCell[,1,drop = F],referencePanel)
sum(test2[,1])

cor(inputCell[,1,drop = F],referencePanel)

test3 <- refCorMerge(mxList,ReferenceList)
test4 <- test3$Cor.merged
View(test4[1:5,1:5])













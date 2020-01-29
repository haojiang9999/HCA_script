## test Function FindCorValueByGenes

x1 <- FindCorValueByGenes(inputCell, referencePanel)
x$singleAbove0$`Fetal__Stomach_5_Pit Progenitor`
x1$singleAbove0_df
colSums(x1$singleAbove0_df)
#### test Function CorValueByGenes_V3

inputCells = GSE97693_Tang_TPM_cells[,1:20] # drop very important
referencePanel = referPanel_of_mix
i=1
x <- CorValueByGenes_V3(inputCells, referencePanel)
x[[1]]
colSums(x[[1]])
table(colSums(x[[1]]) > 0.6)


names(x[[2]])
summary(colSums(x[[2]]))
cor(inputCell, referencePanel)


# verify test
inputCell = GSE97693_Tang_TPM_cells[,2,drop = F] # drop = F to keep rownames
referencePanel = referPanel_of_mix[1:10]
geneList <- intersect(rownames(inputCell),rownames(referencePanel))
inputCell <- inputCell[geneList,,drop=F]
referencePanel <- referencePanel[geneList,]

nonZer0Index <- inputCell > 0                # single cell > 0 genes
#nonZeroInputCell <- inputCell[nonZer0Index,,]
#nonZeroReference <- referencePanel[,i][nonZer0Index]

geneList.2 <- geneList[nonZer0Index]
inputCell[geneList.2,, drop = F]
# calculate corValueGenes



cor(inputCell[geneList.2,, drop = F], referencePanel[geneList.2,])

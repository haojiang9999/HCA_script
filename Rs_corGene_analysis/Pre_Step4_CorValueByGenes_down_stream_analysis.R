######## CorValueByGenes down stream analysis ##########
inputCells = GSE97693_Tang_TPM_cells[,1:20] # drop very important
referencePanel = referPanel_of_mix

CorValueResaults <- CorValueByGenes_V3(inputCells, Filtered_referPanel_of_mix)
#### step1 select high correlation score reference cell types ####
## calculate correlation score of each cell through deferent reference cell types
cellCor <- lapply(CorValueResaults, FUN = colSums)
lapply(cellCor, summary) # summary 
table(cellCor > 0.6)
## corrlation score > X cell tpes
lapply(cellCor, function(x){
  table(x > 0.6)
})
### select top 5 reference cell types for each single cell
lapply(cellCor, function(x){
  head(sort(x, decreasing = T),5)
})
Top5cellTypes <- lapply(cellCor, function(x){
                        names(head(sort(x, decreasing = T),5))
                        })
#### step2 select top 5 reference cell types related genes ####
Top5cellTypes_df <- list()
i = 1
for(i in 1:length(CorValueResaults)){
  Top5cellTypes_df[[i]] <- CorValueResaults[[i]][,Top5cellTypes[[i]]]
}
lapply(Top5cellTypes_df, colSums)  #check

#### step3 find top 200 genes with highest CorValue in top 5 cell types ####
Top5cellTypes_df[[1]]

x <- list()
### generate a table of top genes and top cell types
TopGeneNum = 200
y<- apply(Top5cellTypes_df[[1]],2, function(x){
      TopIndex <- order(x, decreasing = T)[1:TopGeneNum] # how many top genes
      rownames(Top5cellTypes_df[[1]])[TopIndex]
})
x <- Top5cellTypes_df[[1]]
y<- apply(x,2, function(x){
  TopIndex <- order(x, decreasing = T)[1:TopGeneNum] # how many top genes
  rownames(Top5cellTypes_df[[1]])[TopIndex]
})

write.csv(y, file = "test4.csv")

#### genes background for each cell ####
# expressed genes in each cell
# extract background genes from CorValueByGenes_V3
GeneBackGround <- lapply(CorValueResaults, rownames)




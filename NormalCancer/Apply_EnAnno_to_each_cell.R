#### apply EnAnno to each cells
library(parallel)
detectCores()
mc.cores = 6
## Version1: combine all genes from top N cell types.
TopRes[[1]]
as.vector()
test <- mclapply(TopRes, as.vector)
test2 <- mclapply(test, unique)
summary(mclapply(test2, length))
mapply()
lapply(GeneBackGround_ENTREZID, print)
test3 <- cbind(GeneBackGround_ENTREZID, test2)
length(test3)
### convert GeneBackGround_ENTREZID and TopRes to a lapply version
test <- mclapply(TopRes, as.vector)
test2 <- mclapply(test, unique)
ForAnno <- list()
i =1
cellNames <- names(GeneBackGround_ENTREZID)
for(i in 1:length(GeneBackGround_ENTREZID)){
            cell <- list()
            cell$geneMerge <- test2[[i]]
            cell$geneBackGround <- GeneBackGround_ENTREZID[[i]]
            return(cell)
            ForAnno[[cellNames[i]]] <- cell
}

for(i in 1:10){
  print(i)
}
cell <- list()
for(i in 1:10){
  print(i)
  cell$geneMerge <- test2[[i]]
  cell$geneBackGround <- GeneBackGround_ENTREZID[[i]]
  #return(cell)
  ForAnno[[cellNames[i]]] <- cell
  print(i)
  }

mapply(EnAnno, test2, GeneBackGround_ENTREZID)
x <- list(1,1,1,1)
y<- list(2,2,2,2)
mapply(sum, x,y)
assign(GeneBackGround_ENTREZID[[1]])

#### Function prepare for EnAnno
# re-construct data into merged lists
### Inputs were all lists and have same cell order
PreForEnAnno <- function(GeneList, GeneBackGround){
  cell <- list()
  ForAnno <- list()
  cellNames <- names(GeneBackGround)
  for(i in 1:length(GeneBackGround)){
    #print(i)
    cell$geneMerge <- GeneList[[i]]
    cell$geneBackGround <- GeneBackGround_ENTREZID[[i]]
    ForAnno[[cellNames[i]]] <- cell
    print(i)
  }
  return(ForAnno)
}

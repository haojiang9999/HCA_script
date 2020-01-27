#### Function for correlation with reference panel
refCor <- function(input,Reference.list, method = "pearson"){
  corList <- lapply(Reference.list, function(scRef){
    corMatrix( input, scRef)
  })
  corTable <- do.call(rbind.data.frame, corList)
  return(corTable)
}


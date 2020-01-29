##### Function: mean of specific number of columns
sumNcols <- function(df,n){
  y<- list()
  for(i in (0:(dim(df)[2]/n-1)*n+1)){
    
    y[[i]] <- rowMeans(x_sub[,i:(i+n-1)])
    
  }
  Filter(Negate(is.null), y)
}

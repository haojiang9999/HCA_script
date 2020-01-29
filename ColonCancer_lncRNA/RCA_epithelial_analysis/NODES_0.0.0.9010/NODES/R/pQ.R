#' pQ (pseudocounted quantile normalization) for scRNA-seq data.
#' @param data Expression data frame with genes in rows. pQ can be applied to raw-count, FPKM or other expression quantification formats.
#' Gene names should be supplied as rownames whereas sample names as column names. \code{NODES}.
#' @param frac determines a threshold \code{thr = Q1-frac*IQR}; cells with detected genes less than \code{thr+1} would be  thrown away. 
#' @param throw_sd A non zero value will throw away genes that have no non-pseudo count entry after normalization, default is \code{1}.
#' @param hard_outlier A hard lower bound for minimum number of detected genes, default \code{500}.
#' @return normalized data as a  data frame retaining orginal row and column names.
#' @examples 
#' data(data_Trapnell) # load the Trapnell data
#' norm_Data <- pQ(data_Trapnell) # pQ normalization 


pQ <- function(data,frac=0.5,throw_sd=1,hard_outlier=500)
{
  #data = data_Trapnell
  
  # get threshold
  thr = floor(as.numeric(th(data,frac,hard_outlier)))
  
  
  # loading required libraries
  if (!require("preprocessCore")) 
  {
    source("http://bioconductor.org/biocLite.R")
    biocLite("preprocessCore")
    
  }
  require(preprocessCore)
  
  FPKMbackup <- as.matrix(data)
  
  # finding cells which have thr+1 genes
  if(thr != 0)
  {
    bu = apply(FPKMbackup>0,2,function(x) sum(x>0)) >thr
    FPKMreduced<-subset(FPKMbackup,select=bu)
    
    # cells thrown away
    cat(paste("No. of discarded cells:",length(which(bu==FALSE))))
    cat("\n")
    
    
  }else
  {
    FPKMreduced<-FPKMbackup
  }
  

  # minimum number expressed genes in remaining cells  
  det_no <- min(apply(FPKMreduced,2,function(x) sum(x>0)))
  
  # Q normalization
  storage.mode(FPKMreduced) <- "double"
  X <- normalize.quantiles.robust(as.matrix(FPKMreduced),use.median=TRUE,use.log2=FALSE)
  
  rownames(X)<-rownames(FPKMreduced)
  colnames(X)<-colnames(FPKMreduced)
  
  
  # Find value for det_no+1 ranked genes  
  key <- sort(X[,1],decreasing = T)[thr+1]
  
  # level tails
  
  X[X<=key] = key
  
  # throw zero deviation genes
  if(throw_sd!=0)  
  {
    B<-apply(X,1,function(x) sd(x))>0
    X<-X[B,]
    
  }
  
  
  # Report initial and post-processing matrix dimension
  cat(paste("Dimensions of supplied expression matrix:","\n","Genes=",dim(data)[1],";"," Cells=",dim(data)[2],"\n",sep=""))
  cat("\n")
  
  cat(paste("Dimensions of supplied processed matrix:","\n","Genes=",dim(X)[1],";"," Cells=",dim(X)[2],"\n",sep=""))
  cat("\n")
  
  #View(head(X))
  return(X)
  
  
  
}




th <- function(d,frac=0.5,hard_out)
{
  
  # get number of detected genes 
  expr<-apply(d,2,function(x) sum(x>0))
  th = floor(quantile(expr)[2] - abs(quantile(expr)[4]-quantile(expr)[2])*frac)
  
  if (th < hard_out)
  {
    th <- hard_out
  }
  
  th
  
  
}  
  
  
  
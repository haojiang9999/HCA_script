# Function of merged refCor
refCorMerge <- function(mxList, ReferenceList, numCores = 10){
  ### Loading my script
  source("/data8t_4/JH/MyJobs/1_R_script/RCA_seperate_fuction.R")
  source("/data8t_4/JH/MyJobs/1_R_script/refCor.R")
  # ReferenceList =  log2.scReference.list.CV.1500
  #### 1.Correlation of cancer cells with scReference
  ## Using multiple cores
  library(parallel)
  numCores <- numCores
  print(paste("Using",numCores,"parallel cores"))
  Cor.Res<-mclapply(mxList, function(mx){
    refCor(input = mx, Reference.list = ReferenceList, method = "pearson")
  }, mc.cores = numCores)
  Cor.merged <- do.call(cbind, Cor.Res)
  ### recorve the cell names
  CellNames<-mclapply(mxList, function(mx){
    colnames(mx)
  }, mc.cores = numCores)
  CellNames <- as.character(unlist(CellNames))
  colnames(Cor.merged) <- CellNames
  #### 2.Build the sample phenoType infor
  Pheno.merged <- as.data.frame(colnames(Cor.merged))
  colnames(Pheno.merged) <- "Cell_names"
  rownames(Pheno.merged) <- Pheno.merged$Cell_names
  Num<-mclapply(mxList, function(mx){
    ncol(mx)
  }, mc.cores = numCores)
  Num.tb <- do.call(rbind, Num)
  Pheno.merged$Cell_info<- rep(rownames(Num.tb), Num.tb[,1])
  Res <- list(Cor.merged=Cor.merged, Pheno.merged = Pheno.merged)
  return(Res)
  
}



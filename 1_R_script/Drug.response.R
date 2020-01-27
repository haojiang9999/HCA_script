#### Function of dose–response curve to AUC,IC50
Drug.response <- function(DrugRes){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("PharmacoGx")) 
    BiocManager::install("PharmacoGx")
  require(PharmacoGx)
  ## Extract the log10 transformed dose (M)
  drugTable <- DrugRes
  logDose <- sort(drugTable$Unnamed__1)
  logDose <- as.character(logDose)
  drugName <- colnames(drugTable)[1]
  ## Apply calculation to each drugs 
  Drug.res <- apply(drugTable[,-c(1,2)], 2, FUN = function(x){
    #print(x)
    viability <-as.character(rev(x))
    AUC <- computeAUC(logDose, viability,conc_as_log = T)
    IC50 <- computeIC50(logDose, viability, conc_as_log = T)
    Res <- c(AUC,IC50)
    names(Res) <- c(paste0(drugName,"_AUC"),paste0(drugName,"_IC50"))
    return(Res)
  })
  return(Drug.res)
}

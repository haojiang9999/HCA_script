### Calculate the AUC and IC50 using dose–response curve for each organoid
DrugRes

## Extract the log10 transformed dose (M)
drugTable <- DrugRes[[1]]
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

#### Get the AUC,IC50 Resaults
#### Using the function to solve this
source("Drug.response.R")
library(parallel)
#find cores number
detectCores(logical = F)
#set core numbers
mc <- getOption("mc.cores", 10)
x <- mclapply(DrugRes,Drug.response,mc.cores = mc)
x[1:7]
drug.res.P1 <- do.call(rbind.data.frame,x[1:7])
rownames(drug.res.P1) <- unlist(lapply(strsplit(rownames(drug.res.P1),".",fixed = T),"[", 2))
drug.res.P2 <- do.call(rbind.data.frame,x[8:14])
rownames(drug.res.P2) <- unlist(lapply(strsplit(rownames(drug.res.P2),".",fixed = T),"[", 2))
drug.res.P3 <- do.call(rbind.data.frame,x[15:21])
rownames(drug.res.P3) <- unlist(lapply(strsplit(rownames(drug.res.P3),".",fixed = T),"[", 2))
### The resault showed AUC larger the cell more sensitive to the drugs 
### Because the method calculate the 1- viability 
### And the IC50 smaller the cell more sensitive to the drugs
cbind(drug.res.P1,drug.res.P2,drug.res.P3)
drug.res.Description <- "The resault showed AUC larger the cell more sensitive to the drugs, because the method calculate the 1- viability.
And the IC50 smaller the cell more sensitive to the drugs"
drug.res<- cbind(drug.res.P1,drug.res.P2,drug.res.P3)
## only 53 samples had both expression and drug respones data
drug.res <- drug.res[,SamplesBoth]
drug.res<- list(drug.res = drug.res, drug.res.Description = drug.res.Description)


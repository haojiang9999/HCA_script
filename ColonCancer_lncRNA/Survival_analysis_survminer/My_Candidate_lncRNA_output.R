#### My candidate lncRNA Output
Gene_ID <- c("ENSG00000172965","ENSG00000224032",
                 "ENSG00000244479","ENSG00000255198",
                 "ENSG00000225630","ENSG00000245849",
                 "ENSG00000248527","ENSG00000266402")

fileName <- "My_Cand_lncRNA"
#### 3.Data output
######### Survival analysis ##########
dir.create("./Output_survival/")
dir.create(paste0("./Output_survival/",fileName))
setwd(paste0("./Output_survival/",fileName))
### 
library(parallel)
#detectCores()
numCores <- 6
mclapply(Gene_ID,FUN = function(x){
  TCGA.surv(cancerType = "COAD",geneName = x,
            print = F)
}, mc.cores = numCores)
setwd("../")
setwd("../")
######### Expression analysis ##########
dir.create("./Output_expression/")
dir.create(paste0("./Output_expression/",fileName))
setwd(paste0("./Output_expression/",fileName))
mclapply(Gene_ID,FUN = function(x){
  TCGA.expCompare(cancerType = "COAD",geneName = x,
                  print = F)
}, mc.cores = numCores)
### not mc cores work better
lapply(Gene_ID,FUN = function(x){
  TCGA.expCompare(cancerType = "COAD",geneName = x,
                  print = F)
})
setwd("../")
setwd("../")

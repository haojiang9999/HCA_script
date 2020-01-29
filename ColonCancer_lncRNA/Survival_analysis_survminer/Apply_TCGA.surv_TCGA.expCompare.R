############ Gene TCGA survival and expresson analysis
#### 1.Read in data ####
DE_TvN.lnc.DE <- read.csv("dup.lncRNA.csv")
fileName <- "dup"
Gene_ID[2]
Gene_ID <- DE_TvN.lnc.DE$geneName
#### 2.Remove number after '.' 
Gene_ID <- gsub("\\..*","",Gene_ID)
### set up output folder and names


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
setwd("../")
setwd("../")


#### single cell candidate DE lncRNAs' TCGA survival analysis 
epi_stem_DE_TvN.lnc.DE <- read.csv("epi_stem_DE_TvN.lnc.DE.csv")
DE_lncRNAs <- epi_stem_DE_TvN.lnc.DE$geneName
# remove number after '.' 
DE_lncRNAs <- gsub("\\..*","",DE_lncRNAs)
#
table(DE_lncRNAs %in% rownames(TCGA.RSEM.gene.tpm))
### set up output folder
fileName <- "epi_stem"
dir.create(paste0("./Output_survival/",fileName))
setwd(paste0("../Output_survival/",fileName))
### 
library(parallel)
detectCores()
numCores <- 6
mclapply(DE_lncRNAs,FUN = function(x){
  TCGA.surv(cancerType = "COAD",geneName = x,
            print = F)
  }, mc.cores = numCores)
setwd("../")
setwd("../")
dir.create(paste0("./Output_cand_lncRNAs_survival/",fileName))
setwd(paste0("./Output_cand_lncRNAs_survival/",fileName))

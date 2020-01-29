#### Apply_TCGA.expCompare.R ####
#### Step1.Read gene list
#### single cell candidate DE lncRNAs' TCGA survival analysis 
epi_stem_DE_TvN.lnc.DE <- read.csv("epi_stem_DE_TvN.lnc.DE.csv")
DE_lncRNAs <- epi_stem_DE_TvN.lnc.DE$geneName
# remove number after '.' 
DE_lncRNAs <- gsub("\\..*","",DE_lncRNAs)
#
table(DE_lncRNAs %in% rownames(TCGA.RSEM.gene.tpm))
#### Step2.Set up output folder
setwd("./Output_cand_lncRNA_expression/")
#### Step3. Plot figures
library(parallel)
detectCores()
numCores <- 6
mclapply(DE_lncRNAs,FUN = function(x){
  TCGA.expCompare(cancerType = "COAD",geneName = x,
            print = F)
}, mc.cores = numCores)

#### vignettes of NormalCancer ####
### step1 get the cancer cell matrix and reference panel ###
## load reference panel
Filtered_referPanel_of_mix_V1 <- readRDS("test_data/Filtered_referPanel_of_mix_V1.rds")
## load test single cell data
dataSet_for_RCA <- readRDS("test_data/dataSet_for_RCA.rds")
GSE97693_Tang_TPM_cells <- dataSet_for_RCA$GSE97693_Tang_TPM_cells

head(Filtered_referPanel_of_mix_V1)
head(GSE97693_Tang_TPM_cells[1:10,1:10])

### Step2 Calculate CorGene value for each cancer cell in referenece panel ###
CorValueRes <- CorValueByGenes_V3(GSE97693_Tang_TPM_cells, Filtered_referPanel_of_mix_V1)

#CorValueByGenes_V3.688.time <- system.time(CorValueRes <- CorValueByGenes_V3(GSE97693_Tang_TPM_cells, Filtered_referPanel_of_mix_V1))
head(CorValueRes)
CorValueRes[[1]]
### Step3 Top N genes in top N reference cells types
TopRes <- TopCorGenes(CorValueRes, # input list from CorValueByGenes_V3
                      cellTypeNum = 5,  # Number of top correlated cell types
                      TopGeneNum = 200  # Number of top correlated genes in each cell types
                      )
head(TopRes)
### Find the sample names
sampleNames <- names(TopRes$TopCorGenes_mx)
### Step4 Find expression background for each cell
# Genes expressed in all cancer cells
expressIndex <- rowSums(GSE97693_Tang_TPM_cells) > 0
table(expressIndex)
back.list <- names(expressIndex == TRUE)
# convert gene symbol to ENTREZID 
GeneBackGround <- GeneSymbol2GeneID(back.list)
GeneBackGround_ENTREZID <- unique(GeneBackGround$ENTREZID)

### Step5 Gene set enrichmnet analysis
library(parallel)
detectCores()
mc.cores = 6
set.seed( 123, kind = "L'Ecuyer-CMRG" )
####### version1: merge all genes from different normal cell types #########
mergGenes <- mclapply(TopRes, as.vector)
mergGenes <- mclapply(mergGenes, unique)
set.seed(123)
mergGenes_ENTREZID <- mclapply(mergGenes$TopCorGenes_mx, function(x){
                                GeneSymbol2GeneID(x)$ENTREZID
                                },mc.cores = mc.cores)
                                  
mergGenes_ENTREZID <- lapply(mergGenes$TopCorGenes_mx, function(x){
  GeneSymbol2GeneID(x)$ENTREZID
})

names(mergGenes_ENTREZID) <- sampleNames

#### New background version #####
library(parallel)
detectCores()
mc.cores = 25
set.seed( 123, kind = "L'Ecuyer-CMRG" )
time.anno.20cores.688 <- system.time(AnnoRes.688 <- mclapply(mergGenes_ENTREZID, function(x){
                    EnAnno(inputID = x, background = GeneBackGround_ENTREZID, TopN = 10)
                    }, mc.cores = mc.cores,mc.preschedule = FALSE,mc.set.seed = F))
AnnoRes.688[[1]]
### Name the samples
names(AnnoRes.688) <- names(mergGenes_ENTREZID)

#### In here some samples had error happened 
###  need to re-run error ones
### find samples that failed annotation
sample.class<-lapply(AnnoRes.688, class)
fail.samples <- grepl("try",sample.class)
samplesFailed <- names(AnnoRes.688[fail.samples])
library(parallel)
detectCores()
mc.cores = 25
set.seed( 123, kind = "L'Ecuyer-CMRG" )
time.anno.25cores.fialed <- system.time(AnnoRes.688[samplesFailed] <- mclapply(mergGenes_ENTREZID[samplesFailed], function(x){
  EnAnno(inputID = x, background = GeneBackGround_ENTREZID, TopN = 10)
}, mc.cores = mc.cores, mc.preschedule = FALSE,mc.set.seed = F))
## save AnnoRes.688 Normal cancer dataset
AnnoRes.688.merge.datasets <- list(CorValueRes = CorValueRes,
                              TopRes =TopRes,
                              AnnoRes.688.merge = AnnoRes.688)
### Save the dataset 
saveRDS(AnnoRes.688.merge.datasets, file = "AnnoRes.688.merge.datasets.rds")


###############################################################################
#### test version ####
library(parallel)
detectCores()
mc.cores = 20
set.seed( 123, kind = "L'Ecuyer-CMRG" )
time.anno.20cores.688.test <- system.time(AnnoRes.688.test <- mclapply(mergGenes_ENTREZID[1:30], function(x){
  EnAnno(inputID = x, background = GeneBackGround_ENTREZID, TopN = 10)
}, mc.cores = mc.cores, mc.preschedule = FALSE,mc.set.seed = F))

AnnoRes.688.test[[1]]
#### debug ####
debug.test <- lapply(mergGenes_ENTREZID[1:2], function(x){
  EnAnno(inputID = x, background = GeneBackGround_ENTREZID, TopN = 10)
})



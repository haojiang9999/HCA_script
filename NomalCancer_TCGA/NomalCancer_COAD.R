#### TCGA normal cancer ####
# prepare TCGA data
COAD_tpm_symbol

### step1 get the cancer cell matrix and reference panel ###
## load reference panel
Filtered_referPanel_of_mix_V1 <- readRDS("/data8t_4/JH/MyJobs/NormalCancer/test_data/Filtered_referPanel_of_mix_V1.rds")

### Step2 Calculate CorGene value for each cancer cell in referenece panel ###
CorValueRes <- CorValueByGenes_V3(COAD_tpm_symbol, Filtered_referPanel_of_mix_V1)

### Step3 Top N genes in top N reference cells types
TopRes <- TopCorGenes(CorValueRes, # input list from CorValueByGenes_V3
                      cellTypeNum = 5,  # Number of top correlated cell types
                      TopGeneNum = 200  # Number of top correlated genes in each cell types
)
head(TopRes)

### Step4 Find expression background for TCGA COAD
# background genes for COAD
bkGenes <- rownames(COAD_tpm_symbol)[rowSums(COAD_tpm_symbol)>0]
bkGenesID <- GeneSymbol2GeneID(bkGenes)
bkGenesID <- unique(bkGenesID$ENTREZID)

### Step5 gene set enrichment analysis
library(parallel)
detectCores()
mc.cores = 25
set.seed( 123, kind = "L'Ecuyer-CMRG" )
mergGenes <- mclapply(TopRes, as.vector)
mergGenes <- mclapply(mergGenes, unique)
mergGenes_ENTREZID <- mclapply(mergGenes, function(x){
  GeneSymbol2GeneID(input = x, fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Hs.eg.db")$ENTREZID})

time.anno.25cores.COAD <- system.time(AnnoRes.COAD <- mclapply(mergGenes_ENTREZID, function(x){
                                    EnAnno(inputID = x, background = bkGenesID, TopN = 10)
                                      }, mc.cores = mc.cores))

AnnoRes.COAD.flawed <- AnnoRes.COAD
table(lapply(AnnoRes.COAD.flawed, class))
### find samples that failed annotation
x<-lapply(AnnoRes.COAD.flawed, class)
fail.samples <- grepl("try",x)
samplesFailed <- names(AnnoRes.COAD.flawed[fail.samples])
time.anno.25cores.fialed <- system.time(AnnoRes.COAD.flawed <- mclapply(mergGenes_ENTREZID[samplesFailed], function(x){
  EnAnno(inputID = x, background = bkGenesID, TopN = 10)
}, mc.cores = mc.cores))
# re-anno samples by lapply
AnnoRes.COAD.flawed[samplesFailed] <- lapply(mergGenes_ENTREZID[samplesFailed], function(x){
  EnAnno(inputID = x, background = bkGenesID, TopN = 10)})

## save COAD Normal cancer dataset
AnnoRes.COAD.datasets <- list(CorValueRes = CorValueRes,
                              TopRes =TopRes,
                              AnnoRes.COAD.merge = AnnoRes.COAD.flawed)

saveRDS(AnnoRes.COAD.datasets, file = "AnnoRes.COAD.datasets.rds")

###################################################################################


AnnoRes.COAD$TCGA.AZ.4615.01

lapply(mergGenes_ENTREZID$TCGA.AZ.4615.01, function(x){
  EnAnno(inputID = x, background = bkGenesID, TopN = 10)
})


test <- EnAnno(inputID =mergGenes_ENTREZID$TCGA.AZ.4615.01, background = bkGenesID, TopN = 10)

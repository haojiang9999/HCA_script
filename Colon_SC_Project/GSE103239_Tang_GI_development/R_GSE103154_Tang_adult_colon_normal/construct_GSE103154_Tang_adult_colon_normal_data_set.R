# read expression table
filePath <- "/stor/jianghao/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE103154_adult_tissues/GSE103154_All_Merge_umi_tpm_gene.txt.gz"
df.adult.colon.tpm<- read.delim2(gzfile(filePath),row.names = "Gene")
df.adult.colon.tpm[1:10,1:10] # check out
#build sample info
sampleName <- colnames(df.adult.colon.tpm) # extrac sample names
pateintID <-unlist(lapply(strsplit(sampleName,"_"), '[[', 1)) 
tissuePart <- unlist(lapply(strsplit(sampleName,"_"), '[[', 2)) 
tissueType <- unlist(lapply(strsplit(sampleName,"_"), '[[', 3)) 
Batch <- unlist(lapply(strsplit(sampleName,"_"), '[[', 4)) 
cellNumber <- unlist(lapply(strsplit(sampleName,"_"), '[[', 5)) 
pd <- cbind(sampleName, pateintID, tissuePart, tissueType,
            Batch, cellNumber)
rownames(pd) <- sampleName
pd <- as.data.frame(pd)
# feature/gene info
# I do not have it yet
# bulid a dataset
GSE103254.adult.tpm <- list(df.adult.colon.tpm = df.adult.colon.tpm, phenotype = pd)
# file output
saveRDS(GSE103254.adult.tpm, file = "GSE103254.adult.tpm.dataset.rds")





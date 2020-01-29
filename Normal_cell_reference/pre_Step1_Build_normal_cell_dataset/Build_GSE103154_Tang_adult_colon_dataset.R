#### Build_GSE103154_Tang_adult_colon_dataset ####
### Step1 Expression data
filePath <- "/stor/jianghao/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE103154_adult_tissues/GSE103154_All_Merge_umi_tpm_gene.txt.gz"
adult.colon.tpm<- read.table(gzfile(filePath), header = T)
geneNames <- adult.colon.tpm[,1]
rownames(adult.colon.tpm) <- geneNames
adult.colon.tpm <-adult.colon.tpm[,  !(names(adult.colon.tpm) %in% "Gene")]
head(adult.colon.tpm[1:10,1:10])
adult.colon.tpm[1:10,1:10] # check out
summary(colSums(adult.colon.tpm)) 
# the resault showed this was no-transformed tpm value
# But some samples sum was too low not filtered

### Step2 Phenotype information
## extract from sample names
sampleName <- colnames(adult.colon.tpm) # extrac sample names
pateintID <-unlist(lapply(strsplit(sampleName,"_"), '[[', 1)) 
tissuePart <- unlist(lapply(strsplit(sampleName,"_"), '[[', 2)) 
tissueType <- unlist(lapply(strsplit(sampleName,"_"), '[[', 3)) 
Batch <- unlist(lapply(strsplit(sampleName,"_"), '[[', 4)) 
cellNumber <- unlist(lapply(strsplit(sampleName,"_"), '[[', 5)) 
pd <- cbind(sampleName, pateintID, tissuePart, tissueType,
            Batch, cellNumber)
rownames(pd) <- sampleName
pd <- as.data.frame(pd)
## from supplementary data
hcResault <- readRDS("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE103239_Tang_GI_development/R_GSE103154_Tang_adult_colon_normal/Tang_GI_adult_author_cluster_cellType.rds")
### Step3 Gene features
# they don`t provide

### Step4 Build datasets
GSE103154_Tang_adult_colon_datasets <- list(Adult.colon.tpm = adult.colon.tpm,
                                           pd = pd,
                                           hcResault = hcResault)
saveRDS(GSE103154_Tang_adult_colon_datasets, file = "GSE103154_Tang_adult_colon_datasets.rds")





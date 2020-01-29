#### Build GSE95630_Tang_fetal_tissues datasets ####
### Step1 Expression data
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE95630_fetal_tissues/GSE95630_Digestion_TPM_new.txt.gz"
fetal.tissues <- read.table(gzfile(filePath), header = T)
geneNames <- fetal.tissues[,1]
rownames(fetal.tissues) <- geneNames
fetal.tissues <-fetal.tissues[,  !(names(fetal.tissues) %in% "Gene")]
head(fetal.tissues[1:10,1:10])
summary(colSums(fetal.tissues)) 
# the resault showed this was no-transformed tpm value

### Step2 Phenotype information
## extract from sample names
sampleName <- colnames(fetal.tissues)
pd <- as.data.frame(sampleName)
pd$tissue <- unlist(lapply(strsplit(sampleName,"_"), '[[', 1))
pd$emDays <-unlist(lapply(strsplit(sampleName,"_"), '[[', 2))
pd$pateintID <- unlist(lapply(strsplit(sampleName,"_"), '[[', 3))
pd$cellNum <- unlist(lapply(strsplit(sampleName,"_"), '[[', 4))
## from supplementary data
cellTypeCluster <- readRDS("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE103239_Tang_GI_development/GSE95630_fetal_tissues/GSE95630_Tang_fetal_tissues_cellTypeCluster.rds")
### Step3 Gene features
# they don`t provide

### Step4 Build datasets
GSE95630_Tang_fetal_tissues_datasets <- list(Fetal.colon.TPM = fetal.tissues,
                                             pd = pd,
                                             cellTypeCluster = cellTypeCluster)
saveRDS(GSE95630_Tang_fetal_tissues_datasets, 
        file = "GSE95630_Tang_fetal_tissues_datasets.rds")











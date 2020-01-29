filePath <- "/data8t_4/JH/scRNA_seq/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE95630_fetal_tissues/GSE95630_Digestion_TPM_new.txt.gz"
fetal.tissues <- read.table(gzfile(filePath), header = T)
#dim(fetal.tissues)
#head(fetal.tissues)
#rownames(fetal.tissues)
#colnames(fetal.tissues)[1:10]
# add rownames
# phenoData 
geneNames <- fetal.tissues[,1]
rownames(fetal.tissues) <- geneNames
fetal.tissues <-fetal.tissues[,  !(names(fetal.tissues) %in% "Gene")]
sampleName <- colnames(fetal.tissues)
pd <- as.data.frame(sampleName)
pd$tissue <- unlist(lapply(strsplit(sampleName,"_"), '[[', 1))
pd$emDays <-unlist(lapply(strsplit(sampleName,"_"), '[[', 2))
pd$pateintID <- unlist(lapply(strsplit(sampleName,"_"), '[[', 3))
pd$cellNum <- unlist(lapply(strsplit(sampleName,"_"), '[[', 4))
# featureData
geneNames %in% gencode.v19$gene_name
table(geneNames %in% gencode.v19$gene_name)
dim(gencode.v19)
geneNames[!geneNames %in% gencode.v19$gene_name]

# assay
































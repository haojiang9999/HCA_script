#### Build Tang GI development dataset
########################## Adult colon tissue #############################################
## 1.Read expression data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103154
# read expression table
filePath <- "/stor/jianghao/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE103154_adult_tissues/GSE103154_All_Merge_umi_tpm_gene.txt.gz"
df.adult.colon.tpm<- read.delim2(gzfile(filePath),row.names = "Gene")
df.adult.colon.tpm[1:10,1:10] # check out
# convert factor to numeric
df.adult.colon.tpm[] <- lapply(df.adult.colon.tpm, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(df.adult.colon.tpm, class)
summary(colSums(df.adult.colon.tpm))
## 2.Original paper cluster and celltype annotation
# Single-cell multiomics sequencing and analyses of human colorectal cancer
### The authors original cluster resaults
## read file
Adult.hcResault <- read.table("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE103239_Tang_GI_development/R_GSE103154_Tang_adult_colon_normal/Adult_Binary_hc10_result.txt")
Adult.hcResault$hc_10
# check the cluster resaults
library(ggplot2)
p<-ggplot(Adult.hcResault,aes(x=tsne1,y=tsne2,color=as.factor(Adult.hcResault$hc_10)))
p<-p+geom_point(size = 2)
p
# Scatter plot change color
p + scale_color_brewer(palette="Paired")
### 
hcNames <- cbind(seq(10), c("Goblet_1", "Enter_2","Enter_3", "Goblet_2",
                            "Enter_1", "MKI67_High", "Endocrine", "Mesenchymal",
                            "Immune", "OLFM4_High"))
hcNames <- as.data.frame(hcNames)
colnames(hcNames) <- c("hc_10", "cellType")
Adult.hcResault.paper <- Adult.hcResault[,c(1:7,13)]
Adult.hcResault.paper$hc_10 <- as.factor(Adult.hcResault.paper$hc_10)
Adult.hcResault.paper <- dplyr::left_join(Adult.hcResault.paper, hcNames, by = "hc_10")
rownames(Adult.hcResault.paper) <- rownames(Adult.hcResault)
### 3.Gene annotation
Adult.geneNames <- rownames(df.adult.colon.tpm)
source("/data8t_4/JH/MyJobs/1_R_script/GeneSymbol2GeneID.R")
Adult.Gene.Anno <- GeneSymbol2GeneID(as.character(Adult.geneNames), toType = c("ENSEMBL","ENTREZID"))
length(unique(Adult.Gene.Anno$SYMBOL))
## 4.Build Adult colon dataset
Adult_colon_dataset <- list(Adult.colon.exp.TPM = df.adult.colon.tpm, Adult.hcResault.paper = Adult.hcResault.paper,
                            Adult.Gene.Anno = Adult.Gene.Anno)
############################ Fetal tissue ##################################
### 1.1.Read expression data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95630
fetal.GI.tpm <- read.table(gzfile("/data8t_4/JH/scRNA_seq/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE95630_fetal_tissues/GSE95630_Digestion_TPM_new.txt.gz"), header = T,
                           row.names = "Gene")
rownames(fetal.GI.tpm)
colnames(fetal.GI.tpm)[1:10]
summary(colSums(df.adult.colon.tpm))
## 2.Pheno type sample annotation
# add rownames
sampleName <- colnames(fetal.GI.tpm)
Fetal.phenoType <- as.data.frame(sampleName)
Fetal.phenoType$Tissue <- unlist(lapply(strsplit(sampleName,"_"), '[[', 1))
Fetal.phenoType$embryoDays <-unlist(lapply(strsplit(sampleName,"_"), '[[', 2))
Fetal.phenoType$embryoID <- unlist(lapply(strsplit(sampleName,"_"), '[[', 3))
Fetal.phenoType$cellNum <- unlist(lapply(strsplit(sampleName,"_"), '[[', 4))
## 3.Original paper cluster Annotation
Esophagus_KNN <- read.csv("Esophagus_KNN_Result.csv")
LIntes_KNN <- read.csv("L-Intes_KNN_Result.csv")
SIntes_KNN <- read.csv("S-Intes_KNN_Result.csv")
Stomach_KNN <- read.csv("Stomach_KNN_Result.csv")
Feathers_organs <- read.csv("Table_1_Features_of_fetal_digestive_organs.csv")
Fetal.Cluster.paper <- list(Esophagus_KNN = Esophagus_KNN,
                            LIntes_KNN = LIntes_KNN,
                            SIntes_KNN = SIntes_KNN,
                            Stomach_KNN = Stomach_KNN,
                            Feathers_organs = Feathers_organs)

## 4.Gene annotation 
Fetal.geneNames <- rownames(fetal.GI.tpm)
source("/data8t_4/JH/MyJobs/1_R_script/GeneSymbol2GeneID.R")
Fetal.Gene.Anno <- GeneSymbol2GeneID(as.character(Fetal.geneNames), toType = c("ENSEMBL","ENTREZID"))
length(unique(Fetal.Gene.Anno$SYMBOL))
## 5.Build Fetal GI develepment dataset
Fetal_GI_dataset <- list(Fetal.GI.exp.TPM = fetal.GI.tpm, Fetal.phenoType = Fetal.phenoType,
                         Fetal.Cluster.paper = Fetal.Cluster.paper, Fetal.Gene.Anno = Fetal.Gene.Anno)



#### Save the data set ####
saveRDS(Adult_colon_dataset, file = "Tang_Adult_colon_dataset.rds")
saveRDS(Fetal_GI_dataset, file = "Tang_Fetal_GI_dataset.rds")

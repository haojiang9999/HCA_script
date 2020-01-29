#### Step1 Read data from different paper ####
### read csv files from 
# Automated brightfeld morphometry of 3D organoid populations by OrganoSeg
## All samples they sequenced
MOESM4 <- read.csv("41598_2017_18815_MOESM4_ESM.csv")
## Organoid samples filter gene names and re-normalized TPM (subset of MOESM4)
MOESM5 <- read.csv("41598_2017_18815_MOESM5_ESM.csv")
summary(colSums(MOESM5[,-(1:2)]))
## Check the sample names
colnames(MOESM5)
colnames(MOESM4)
table(colnames(MOESM5) %in% colnames(MOESM4))
GeneAnno <- MOESM5[,c(1,2)]
GeneAnno_Full <- MOESM4[,(1:8)]
expMatrix <- MOESM5[,-(1:2)]
rownames(expMatrix) <- GeneAnno[,1]
colnames()
# read files from paper2
# Intra-tumour diversification in colorectal cancer at the single-cell level
csvs <- list.files("MOESM10_ESM_CSV/")
csvs <- paste0("/data8t_4/JH/MyJobs/Read_dataset/3D_organoid/MOESM10_ESM_CSV/", csvs)
CSVfiles = lapply(csvs, read.csv)
#CSVfiles[[2]]
#x<-CSVfiles[[9]]
#### Step2 build the drug response table ####
## Calculate the average survival ratio
DrugRes<- lapply(CSVfiles, function(x){
  colnames(x) <- gsub(".","_",colnames(x), fixed = T)
  # collect colnames
  sampleNames <- colnames(x[,-(1:2)])[seq(1, ncol(x[,-(1:2)]), 4)]
  drug <- x[,(1:2)]
  x_sub <<- x[,-(1:2)]
  meanSur <- meanNcols(df = x_sub, n = 4)
  mx<- do.call("cbind", meanSur)
  colnames(mx) <- sampleNames
  
  print("1")
  #return(df1)
  cbind(drug,mx)
})

#DrugRes[[2]]
### Add names for DrugRes
SampleNames<- lapply(DrugRes, function(x){
  sampleNames <- colnames(x)
  sName <- paste0(sampleNames[1],"_",sampleNames[3])
  return(sName)
})
names(DrugRes) <- SampleNames
### Find the sample have expression data and drug response data
# Every sample organoid had 7 drug response data
grep("P2_T1_3",DrugRes)
sampleNames <- unlist(lapply(CSVfiles, colnames))
DrugSampleNames <- sampleNames[grep("P",sampleNames)]
DrugSampleNames <- gsub(".","_",DrugSampleNames, fixed = T)
# only 53 samples had both expression and drug respones data
table(colnames(MOESM5) %in% DrugSampleNames)
SamplesBoth <- colnames(MOESM5)[colnames(MOESM5) %in% DrugSampleNames]
expMatrix_both <- expMatrix[,SamplesBoth]
summary(colSums(expMatrix_both))
### build colon Organoid dataset
Colon_Orgnoid_dataset<- list(Orgnoid_expression.TPM = expMatrix_both,
                             Drug_Respones = DrugRes,
                             GeneAnno = GeneAnno,
                             GeneAnno_Full = GeneAnno_Full,
                             Drug_AUC_IC50 = drug.res # from PharmacoGx
                             )
saveRDS(Colon_Orgnoid_dataset, file = "Colon_Orgnoid_dataset.rds")


################ Convert exp Matrix ENSEMBL ID to SYMBOL ###################
#expMatrix_both[1:5,1:5]
#ENSEMBL <- rownames(expMatrix_both)
#source("/data8t_4/JH/MyJobs/1_R_script/EnsemblID2GeneSymbol.R")
#eg <- EnsemblID2GeneSymbol(ENSEMBL)
## remove EnsemblID not have Symbol
#eg.sub <- eg[nchar(eg$hgnc_symbol) != 0,]
#rownames(eg.sub) <- eg.sub$ensembl_gene_id
## 1.select ENSEMBL had a SYMBOL paire
#ENSEMBL.dup <- eg.sub[duplicated(eg.sub$hgnc_symbol),]$ensembl_gene_id
### Step1.Build Non-dup symbol expression matrix
# Symbol of dup or non-dup genes
#Dup.genes <- as.character(unique(eg.sub[duplicated(eg.sub$hgnc_symbol),]$hgnc_symbol))
#NonDup.genes <-as.character(unique(GeneAnno.sub[!duplicated(GeneAnno.sub$GeneSymbol),]$GeneSymbol))
#NonDup.genes <-as.character(eg.sub[!eg.sub$hgnc_symbol %in% Dup.genes,]$hgnc_symbol)
### Expression dataframe of Non-dup genes
#Index.NonDup<- eg.sub[!eg.sub$hgnc_symbol %in% Dup.genes,]$ensembl_gene_id## remove dup genes
#table(duplicated(Index.NonDup))
#expMatrix_both.uni1 <- expMatrix_both[Index.NonDup,]
## change ensembl_gene_id to hgnc_symbol
#rownames(expMatrix_both.uni1) <- eg.sub[Index.NonDup,]$hgnc_symbol
### Step2.Build duplicated ensembl_gene_id expression matrix using Max expression one
#expMatrix_both.uni2 <- sapply(Dup.genes, FUN = function(symbol){
  #symbol= "VPS52"
  Index.Dup<- rownames(eg.sub[eg.sub$hgnc_symbol %in% symbol,])
  exp.uni2 <- expMatrix_both[Index.Dup,]
  # Find Ensembl ID that have highest sum value
  exp.rowsum <- rowSums(exp.uni2)
  Name.Max <- names(exp.rowsum[exp.rowsum == max(exp.rowsum)])
  exp.uni2 <- exp.uni2[Name.Max,]
})
#expMatrix_both.uni2 <- t(expMatrix_both.uni2)
#### Step3.Combine two unique table
#expMatrix_both.uniGene <- rbind(expMatrix_both.uni1,expMatrix_both.uni2)
#sapply(expMatrix_both.uniGene, class)
# convert list to numeric
#expMatrix_both.uniGene <- sapply(expMatrix_both.uniGene, unlist)
#expMatrix_both.uniGene<- as.data.frame(expMatrix_both.uniGene)
#### Convert exp Matrix ENSEMBL ID to SYMBOL
source("/data8t_4/JH/MyJobs/1_R_script/ExpMxID2Symbol.R")
expMatrix_both.uniGene <- ExpMxID2Symbol(expMatrix_both)
sapply(expMatrix_both.uniGene, class)
rownames(expMatrix_both.uniGene)
Colon_Orgnoid_uniGene_dataset<- list(Orgnoid_expression.uniGene.TPM = expMatrix_both.uniGene,
                             Drug_Respones = DrugRes,
                             GeneAnno = GeneAnno,
                             GeneAnno_Full = GeneAnno_Full,
                             Drug_AUC_IC50 = drug.res # from PharmacoGx
)
saveRDS(Colon_Orgnoid_uniGene_dataset, file = "Colon_Orgnoid_uniGene_dataset.rds")










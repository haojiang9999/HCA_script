### 4.Read_COAD_TCGA_mutation_maf_files.R
# Download from TCGA data portal
# https://portal.gdc.cancer.gov/files/70cb1255-ec99-4c08-b482-415f8375be3f
# https://portal.gdc.cancer.gov/files/8177ce4f-02d8-4d75-a0d6-1c5450ee08b0
# https://portal.gdc.cancer.gov/files/faa5f62a-2731-4867-a264-0e85b7074e87
## Clean header with # 
## and gunzip file can read  OK
### Step1.read data using maftools
library(maftools)
TCGA.COAD.mutect <- "/stor/jianghao/TCGA/COAD/maf/TCGA.COAD.mutect.somatic.maf"
COAD.mutect <- read.maf(maf = TCGA.COAD.mutect)
#Shows sample summry.
getSampleSummary(COAD.mutect)
#Shows all fields in MAF
getFields(COAD.mutect)
 oncostrip(maf = COAD.mutect, genes = c('KRAS','BRAF', 'APC','TP53'))

oncoplot(maf = COAD.mutect, top = 10)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = COAD.mutect, basename = 'TCGA.COAD.mutect')
### Step2.Read summary files
TCGA.COAD.mutect_geneSummary <- read.table("TCGA.COAD.mutect_geneSummary.txt",header = T)
TCGA.COAD.mutect_sampleSummary <- read.table("TCGA.COAD.mutect_sampleSummary.txt",header = T)
TCGA.COAD.mutect_summary <- read.table("TCGA.COAD.mutect_summary.txt",header = T)
### change TCGA.COAD.mutect_sampleSummary rownames
sampleID <- as.character(TCGA.COAD.mutect_sampleSummary$Tumor_Sample_Barcode)
sampleID <- substr(sampleID,1,nchar(sampleID)-13)
rownames(TCGA.COAD.mutect_sampleSummary) <- gsub("-",".",sampleID)

### Step3.Mutation genes by samples
TCGA.COAD.mutect_mutCountMatrix <- mutCountMatrix(COAD.mutect, removeNonMutated = F)
## Change sample names
sampleID <- colnames(TCGA.COAD.mutect_mutCountMatrix)
sampleID <- substr(sampleID,1,nchar(sampleID)-13)
colnames(TCGA.COAD.mutect_mutCountMatrix) <- gsub("-",".",sampleID)

#### Step4.Extract Information I needed
#1)Total mutation
tatolMut.df <- data.frame(rownames = rownames(TCGA.COAD.mutect_sampleSummary),
                          totalMut = TCGA.COAD.mutect_sampleSummary$total)
#2)Genes whether mutated in samples?
Genes <- c('KRAS','BRAF', 'APC','TP53')
GeneMut.df <- TCGA.COAD.mutect_mutCountMatrix[Genes,]
### Convert mutate to 1 unmutate to 0
GeneMut.df[GeneMut.df>0]<-1 
GeneMut.df <- t(GeneMut.df)
GeneMut.df <- data.frame(rownames = rownames(GeneMut.df), GeneMut.df)
#3)Merge two dataframe
TCGA.COAD.Mut <- dplyr::left_join(GeneMut.df, tatolMut.df, by = "rownames")
TCGA.COAD.Mut$rownames
TCGA.COAD.Mut$patient_barcode <- substr(TCGA.COAD.Mut$rownames,1,nchar(as.character(TCGA.COAD.Mut$rownames))-3)




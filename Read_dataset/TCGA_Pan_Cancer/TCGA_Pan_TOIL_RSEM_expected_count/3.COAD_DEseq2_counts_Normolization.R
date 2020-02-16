#### 3.COAD_DEseq2_counts_Normolization.R
### 1.Read data
cancerType = "COAD"
filePath <- paste0("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_Pan_TOIL_RSEM_expected_count/",cancerType,"_UCSC_Toil_expected_count_dataset.rds")
TCGA.Toil.expected_count.dataset <- readRDS(filePath)
TCGA.RSEM.gene.expected_count <- TCGA.Toil.expected_count.dataset[[paste0(cancerType,".RSEM.gene.expected_count_round")]]
TCGA.pheno <- TCGA.Toil.expected_count.dataset[[paste0(cancerType,".pheno")]]
TCGA.RSEM.gene.expected_count[1:5,1:5]
dim(TCGA.RSEM.gene.expected_count)
#### 2. Select tumor samples
## sample type ###
sampleName <- rownames(TCGA.pheno)
samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
table(samplType)
## "01" was primary tumor, "11' was normal tissue
tumorIndex <- samplType == "01"
#normalIndex <- samplType == "11"
TCGA.pheno.tumor <- TCGA.pheno[tumorIndex,]
tumorSamples <- rownames(TCGA.pheno.tumor)
TCGA.RSEM.gene.expected_count.Tumor <- TCGA.RSEM.gene.expected_count[,tumorSamples]
TCGA.RSEM.gene.expected_count.Tumor[1:5,1:5]
#### 3.DE gene analysis using DEseq2
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(DESeq2)
library("BiocParallel")
register(MulticoreParam(5))
dds.tumor <- DESeqDataSetFromMatrix(countData = TCGA.RSEM.gene.expected_count.Tumor,
                              colData = TCGA.pheno.tumor[colnames(TCGA.RSEM.gene.expected_count.Tumor),],)

























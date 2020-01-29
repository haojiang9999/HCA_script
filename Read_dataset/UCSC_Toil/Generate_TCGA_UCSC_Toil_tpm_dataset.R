#### 1.read RNA-seq expression data sets ####
filePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/gene_expression_RNAseq/tcga_RSEM_gene_tpm/tcga_RSEM_gene_tpm.gz"
#### Data (file names: *.rsem_genes.results) are downloaded, 
#    tpm values are extracted, log2(x+0.001) transformed
tcga_RSEM_gene_tpm <- read.table(filePath, header = T)
## add first column to rownames
rownames(tcga_RSEM_gene_tpm) <- tcga_RSEM_gene_tpm$sample
# remove first columns
tcga_RSEM_gene_tpm <- tcga_RSEM_gene_tpm[,-1]
head(tcga_RSEM_gene_tpm[1:10,1:10])
saveRDS(tcga_RSEM_gene_tpm, file = "tcga_RSEM_gene_tpm.rds")
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read gencode.v23.annotation.gene.probemap ####
geneAnnoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/gene_expression_RNAseq/tcga_RSEM_gene_tpm/gencode.v23.annotation.gene.probemap"
gencode.v23.annotation <- read.table(geneAnnoFilePath)


#### 4.Generate TCGA_UCSC_Toil_tpm_dataset ####
#i="BRCA"

for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- phenotype_sp$cancer.type.abbreviation == i
  TCGA.pheno <- phenotype_sp[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$sample)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  rownames(TCGA.pheno) <- TCGA.sampleID
  ## Step3 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(tcga_RSEM_gene_tpm) %in% TCGA.sampleID
  TCGA_RSEM_gene_tpm  <- tcga_RSEM_gene_tpm[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene_tpm),]
  ## Step4 Reverse log2 + 0.001 expression
  TCGA_RSEM_gene_tpm<-apply(TCGA_RSEM_gene_tpm,2, function(x){
    2^x - 0.001
  })
  # replace negative value to zeros
  TCGA_RSEM_gene_tpm[TCGA_RSEM_gene_tpm <0 ] <- 0 
  ## Step5 Biuld TCGA data sets
  TCGA_UCSC_Toil_tpm <- list(exp = TCGA_RSEM_gene_tpm,
                             pheno = TCGA.pheno.exp,
                             gencode.v23.annotation = gencode.v23.annotation)
  names(TCGA_UCSC_Toil_tpm)<-c(paste0(i,".RSEM.gene.tpm"),pheno = paste0(i,".pheno"), "gencode.v23.annotation")
  saveRDS(TCGA_UCSC_Toil_tpm, file = paste0(i,"_UCSC_Toil_tpm_dataset.rds"))
}




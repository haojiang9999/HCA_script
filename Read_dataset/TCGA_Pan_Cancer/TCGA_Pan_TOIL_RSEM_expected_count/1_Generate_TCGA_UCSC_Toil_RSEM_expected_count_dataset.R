#### 1.read RNA-seq expression data sets ####
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_expected_count/tcga_gene_expected_count.gz"
#### https://xenabrowser.net/datapages/?dataset=tcga_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Data (file names: *.rsem_genes.results) are downloaded, expected_count values are extracted, log2(x+1) transformed, and combined.
tcga_RSEM_expected_count <- read.table(filePath, header = T)
rownames(tcga_RSEM_expected_count) <- tcga_RSEM_expected_count$sample
# remove first columns
tcga_RSEM_expected_count <- tcga_RSEM_expected_count[,-1]
head(tcga_RSEM_expected_count[1:10,1:10])
saveRDS(tcga_RSEM_expected_count, file = "tcga_RSEM_expected_count.rds")
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read gencode.v23.annotation.gene.probemap ####
geneAnnoFilePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_expected_count/gencode.v23.annotation.gene.probemap"
gencode.v23.annotation <- read.table(geneAnnoFilePath)
#### 4.Generate TCGA_UCSC_Toil_tpm_dataset ####
#i="ACC"

for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- phenotype_sp$cancer.type.abbreviation == i
  TCGA.pheno <- phenotype_sp[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$sample)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  rownames(TCGA.pheno) <- TCGA.sampleID
  ## Step3 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(tcga_RSEM_expected_count) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- tcga_RSEM_expected_count[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Reverse log2 (x+ 1) expression
  TCGA_RSEM_gene<-apply(TCGA_RSEM_gene,2, function(x){
    2^x - 1
  })
  TCGA_RSEM_gene <- round(TCGA_RSEM_gene)
  # replace negative value to zeros
  TCGA_RSEM_gene[TCGA_RSEM_gene <0 ] <- 0 
  ## Step5 Biuld TCGA data sets
  TCGA_UCSC_Toil_tpm <- list(exp = TCGA_RSEM_gene,
                             pheno = TCGA.pheno.exp,
                             gencode.v23.annotation = gencode.v23.annotation)
  names(TCGA_UCSC_Toil_tpm)<-c(paste0(i,".RSEM.gene.expected_count_round"),pheno = paste0(i,".pheno"), "gencode.v23.annotation")
  saveRDS(TCGA_UCSC_Toil_tpm, file = paste0(i,"_UCSC_Toil_expected_count_dataset.rds"))
}


#### 1.read RNA-seq expression data sets ####
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_norm_count/tcga_RSEM_Hugo_norm_count.gz"
#### https://xenabrowser.net/datapages/?dataset=tcga_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Data (file names: *.rsem.genes.norm_counts.hugo.tab) are log2(x+1) transformed, and combined
tcga_RSEM_norm_count <- read.table(filePath, header = T)
## add first column to rownames
rownames(tcga_RSEM_norm_count) <- tcga_RSEM_norm_count$sample
# remove first columns
tcga_RSEM_norm_count <- tcga_RSEM_norm_count[,-1]
head(tcga_RSEM_norm_count[1:10,1:10])
saveRDS(tcga_RSEM_norm_count, file = "tcga_RSEM_norm_count.rds")
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read gencode.v23.annotation.gene.probemap ####
geneAnnoFilePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_norm_count/hugo_gencode_good_hg38_v23comp_probemap"
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
  TCGA.sampleID.exp <- colnames(tcga_RSEM_norm_count) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- tcga_RSEM_norm_count[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Reverse log2(x+1) expression
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
  names(TCGA_UCSC_Toil_tpm)<-c(paste0(i,".RSEM.gene.norm_count_round"),pheno = paste0(i,".pheno"), "gencode.v23.annotation")
  saveRDS(TCGA_UCSC_Toil_tpm, file = paste0(i,"_UCSC_Toil_norm_count_dataset.rds"))
}
################ Log2(x+1) #######################
#log2(norm_count+1)
# https://xenabrowser.net/datapages/?dataset=tcga_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Data (file names: *.rsem.genes.norm_counts.hugo.tab) are log2(x+1) transformed, and combined

for(i in TCGA.cancer.types){
  ### Step1 separate data by cancer types ###
  TCGA_Index <- phenotype_sp$cancer.type.abbreviation == i
  TCGA.pheno <- phenotype_sp[TCGA_Index,]
  #### Step2 Convert and add rownames to pheno table #### 
  TCGA.sampleID <- as.character(TCGA.pheno$sample)
  TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
  rownames(TCGA.pheno) <- TCGA.sampleID
  ## Step3 Extract expression data by pheno data
  TCGA.sampleID.exp <- colnames(tcga_RSEM_norm_count) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- tcga_RSEM_norm_count[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Reverse log2(x+1) expression
  #TCGA_RSEM_gene<-apply(TCGA_RSEM_gene,2, function(x){
  #  2^x - 1
  #})
  #TCGA_RSEM_gene <- round(TCGA_RSEM_gene)
  # replace negative value to zeros
  #TCGA_RSEM_gene[TCGA_RSEM_gene <0 ] <- 0 
  ## Step5 Biuld TCGA data sets
  TCGA_UCSC_Toil_RSEM <- list(exp = TCGA_RSEM_gene,
                             pheno = TCGA.pheno.exp,
                             gencode.v23.annotation = gencode.v23.annotation)
  names(TCGA_UCSC_Toil_RSEM)<-c(paste0(i,".RSEM.gene.norm_count_Log2(x+1)"),pheno = paste0(i,".pheno"), "gencode.v23.annotation")
  saveRDS(TCGA_UCSC_Toil_RSEM, file = paste0(i,"_UCSC_RSEM_norm_count_Log2(x+1)_hugo_dataset.rds"))
}

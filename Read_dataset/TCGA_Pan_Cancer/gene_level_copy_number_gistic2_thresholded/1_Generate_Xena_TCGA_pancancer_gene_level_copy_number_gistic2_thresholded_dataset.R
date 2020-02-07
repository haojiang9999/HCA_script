#### 1.Generate_Xena_TCGA_pancancer_gene_level_copy_number_gistic2_thresholded_dataset.R
#### 1.read gene_level_copy_number_gistic2_thresholded_n_10845 ####
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/Copy_number/gene_level_copy_number_gistic2_thresholded_n_10845/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
#### https://xenabrowser.net/datapages/?dataset=TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
gene_level_copy_number_gistic2_thresholded <- read.table(filePath, header = T)
gene_level_copy_number_gistic2_thresholded[1:5,1:5]
## add first column to rownames
rownames(gene_level_copy_number_gistic2_thresholded) <- gene_level_copy_number_gistic2_thresholded$Sample
# remove first columns
gene_level_copy_number_gistic2_thresholded <- gene_level_copy_number_gistic2_thresholded[,-1]
head(gene_level_copy_number_gistic2_thresholded[1:10,1:10])
saveRDS(gene_level_copy_number_gistic2_thresholded, file = "gene_level_copy_number_gistic2_thresholded.rds")
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read gencode.v23.annotation.gene.probemap ####
geneAnnoFilePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_norm_count/hugo_gencode_good_hg38_v23comp_probemap"
gencode.v23.annotation <- read.table(geneAnnoFilePath)
#### 4.Generate Xena_TCGA_copynumber_gene_gistic2_thresholded_dataset ####
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
  TCGA.sampleID.exp <- colnames(gene_level_copy_number_gistic2_thresholded) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- gene_level_copy_number_gistic2_thresholded[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Biuld TCGA data sets
  TCGA_UCSC_gene_gistic2_thresholded <- list(gene_copynumber_gistic2_thresholded = TCGA_RSEM_gene,
                             pheno = TCGA.pheno.exp,
                             gencode.v23.annotation = gencode.v23.annotation)
  names(TCGA_UCSC_gene_gistic2_thresholded)<-c(paste0(i,".gene.copynumber.gistic2.thresholded"),pheno = paste0(i,".pheno"), "gencode.v23.annotation")
  saveRDS(TCGA_UCSC_gene_gistic2_thresholded, file = paste0(i,"_UCSC_gene_gistic2_thresholded_dataset.rds"))
}









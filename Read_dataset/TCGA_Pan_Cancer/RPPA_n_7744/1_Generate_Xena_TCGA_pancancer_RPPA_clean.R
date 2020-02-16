#### 1_Generate_Xena_TCGA_pancancer_RPPA_clean.R
#### 1.read Xena_TCGA_pancancer_RPPA_clean_n_7744 ####
#### https://xenabrowser.net/datapages/?dataset=TCGA-RPPA-pancan-clean.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/RPPA_n_7744/TCGA-RPPA-pancan-clean.xena.gz"
TCGA_RPPA_pancan_clean.xena <- read.table(filePath ,header = T,row.names = "SampleID")
TCGA_RPPA_pancan_clean.xena[1:5,1:5]
saveRDS(TCGA_RPPA_pancan_clean.xena, file = "TCGA_RPPA_pancan_clean.xena.rds")
#### 2.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
## How many cancer types
TCGA.cancer.types <- names(table(phenotype_sp$cancer.type.abbreviation))
#### 3. Read gencode.v23.annotation.gene.probemap ####
geneAnnoFilePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_norm_count/hugo_gencode_good_hg38_v23comp_probemap"
gencode.v23.annotation <- read.table(geneAnnoFilePath)
#### 4.Generate_Xena_TCGA_pancancer_RPPA_clean_dataset ####
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
  TCGA.sampleID.exp <- colnames(TCGA_RPPA_pancan_clean.xena) %in% TCGA.sampleID
  TCGA_RSEM_gene  <- TCGA_RPPA_pancan_clean.xena[,TCGA.sampleID.exp]
  TCGA.pheno.exp <- TCGA.pheno[colnames(TCGA_RSEM_gene),]
  ## Step4 Biuld TCGA data sets
  TCGA_RPPA_clean.xena <- list(TCGA_RPPA_pancan_clean.xena = TCGA_RSEM_gene,
                                             pheno = TCGA.pheno.exp,
                                             gencode.v23.annotation = gencode.v23.annotation)
  names(TCGA_RPPA_clean.xena)<-c(paste0(i,".RPPA.pancan.clean.xena"),pheno = paste0(i,".pheno"), "gencode.v23.annotation")
  saveRDS(TCGA_RPPA_clean.xena, file = paste0(i,"_TCGA_RPPA_clean_xena_dataset.rds"))
}


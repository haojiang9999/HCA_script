#### 1.Generation_COAD_RNAseqV2_dataset.R
# https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Level_3 data (file names: *.rsem.genes.normalized_results) are downloaded from TCGA DCC, 
# log2(x+1) transformed,and processed at UCSC into Xena repository
### 1.Read expression data 20530 * 329
COAD.HiSeqV2 <- read.table("/stor/jianghao/Xena/TCGA_Hub/TCGA_Colon_Cancer_COAD/gene_expression_RNAseq/IlluminaHiSeq (n=329) TCGA Hub/HiSeqV2.gz",
                           header = T, row.names = "sample")
COAD.HiSeqV2[1:5,1:5]
class(COAD.HiSeqV2$TCGA.CA.5256.01)
### 1.1 Read metadata
library(jsonlite)
COAD.HiSeqV2.metadata <- fromJSON("/stor/jianghao/Xena/TCGA_Hub/TCGA_Colon_Cancer_COAD/gene_expression_RNAseq/IlluminaHiSeq (n=329) TCGA Hub/HiSeqV2.json")
### 2.Gene annotation data
hugo_gencode_good_hg19_V24lift37_probemap <- read.table("/stor/jianghao/Xena/TCGA_Hub/TCGA_Colon_Cancer_COAD/gene_expression_RNAseq/IlluminaHiSeq (n=329) TCGA Hub/hugo_gencode_good_hg19_V24lift37_probemap",
                                                        header = T)
hugo_gencode_good_hg19_V24lift37_probemap[1:5,1:5]
### 2.1 read annotation metadata
hugo_gencode_good_hg19_V24lift37.metadata <- fromJSON("/stor/jianghao/Xena/TCGA_Hub/TCGA_Colon_Cancer_COAD/gene_expression_RNAseq/IlluminaHiSeq (n=329) TCGA Hub/hugo_gencode_good_hg19_V24lift37_probemap.json")

### 3.Read phenotype data #### 
phenoFilePath <- "/stor/jianghao/Xena/UCSC_Toil/TCGA_Pan_Cancer_PANCAN/phenotype/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
phenotype_sp <- read.delim(phenoFilePath) 
table(phenotype_sp$cancer.type.abbreviation)
phenotype_sp$sample
TCGA.sampleID <- as.character(phenotype_sp$sample)
TCGA.sampleID <- gsub("-",".",TCGA.sampleID)
rownames(phenotype_sp) <- TCGA.sampleID
# subset COAD data
COAD.pheno <- phenotype_sp[colnames(COAD.HiSeqV2),]

#### 4.Read gencode.v23.annotation.gene.probemap ####
geneAnnoFilePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/gene_expression_RNAseq/TOIL_RSEM_expected_count/gencode.v23.annotation.gene.probemap"
gencode.v23.annotation <- read.table(geneAnnoFilePath,header = T)

### 5.Save dataset
COAD_RNAseqV2_dataset_log2 <- list(COAD.HiSeqV2.log2 =COAD.HiSeqV2,
                                   COAD.HiSeqV2.metadata = COAD.HiSeqV2.metadata,
                                  hugo_gencode_good_hg19_V24lift37_probemap=hugo_gencode_good_hg19_V24lift37_probemap,
                                  hugo_gencode_good_hg19_V24lift37.metadata = hugo_gencode_good_hg19_V24lift37.metadata,
                                  COAD.pheno=COAD.pheno,
                                  gencode.v23.annotation=gencode.v23.annotation)
saveRDS(COAD_RNAseqV2_dataset_log2, file = "TCGA_COAD_RNAseqV2_normalized_log2_dataset.rds")




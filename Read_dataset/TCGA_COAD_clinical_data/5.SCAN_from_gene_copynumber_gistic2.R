#### 5.SCAN_from_gene_copynumber_gistic2.R
## COAD_UCSC_gene_gistic2_thresholded_dataset
# Somatic copy number alterations (SCNA) 
### 1.read the COAD_UCSC_gene_gistic2_thresholded_dataset from UCSC XENA
COAD_gene_gistic2_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/gene_level_copy_number_gistic2_thresholded/COAD_UCSC_gene_gistic2_thresholded_dataset.rds")
COAD.gene.copynumber.gistic2.thresholded <- COAD_gene_gistic2_dataset$COAD.gene.copynumber.gistic2.thresholded
### 2.Counts the SCNA like paper
# The consensus molecular subtypes of colorectal cancer doi:10.1038/nm.3967
# We counted GISTIC scores −2/−1/+1/+2 as events for SCNA estimation 
#(<Q1 was considered low and ≥Q1 was considered high).
COAD.SCNA.gistic2 <- COAD.gene.copynumber.gistic2.thresholded
COAD.SCNA.gistic2[COAD.SCNA.gistic2 !=0 ] <- 1 # all non-zero were copy number alterations
### 3.Count SCNA for each samples
COAD.SCNA.sample.counts <- t(COAD.SCNA.gistic2)
COAD.SCNA.sample.counts <- rowSums(COAD.SCNA.sample.counts)
COAD.SCNA.sample.counts <- data.frame(rownames = names(COAD.SCNA.sample.counts),
                                      SCNA.gene.counts = COAD.SCNA.sample.counts)







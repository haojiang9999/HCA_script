#### 2.Read_rds_test.R
ACC <- readRDS("ACC_UCSC_gene_gistic2_thresholded_dataset.rds")
ACC.gene.copynumber.gistic2.thresholded <- ACC$ACC.gene.copynumber.gistic2.thresholded
ACC.pheno <- ACC$ACC.pheno
gencode.v23.annotation <- ACC$gencode.v23.annotation
ACC$gistic2.thresholded.metadata

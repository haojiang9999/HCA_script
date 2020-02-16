#### 2.test
ACC <- readRDS("ACC_TCGA_RPPA_clean_xena_dataset.rds")
ACC.RPPA.pancan.clean.xena <- ACC$ACC.RPPA.pancan.clean.xena
ACC.pheno <- ACC$ACC.pheno
gencode.v23.annotation <- ACC$gencode.v23.annotation

COAD <- readRDS("COAD_TCGA_RPPA_clean_xena_dataset.rds")
COAD.RPPA.pancan.clean.xena  <- COAD$COAD.RPPA.pancan.clean.xena

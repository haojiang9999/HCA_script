#### 2.test
ACC <- readRDS("ACC_TCGA_RPPA_clean_xena_dataset.rds")
ACC.RPPA.pancan.clean.xena <- ACC$ACC.RPPA.pancan.clean.xena
ACC.pheno <- ACC$ACC.pheno
gencode.v23.annotation <- ACC$gencode.v23.annotation

COAD <- readRDS("COAD_TCGA_RPPA_clean_xena_dataset.rds")
COAD.RPPA.pancan.clean.xena  <- COAD$COAD.RPPA.pancan.clean.xena
COAD.RPPA.pancan.clean.xena
dim(COAD.RPPA.pancan.clean.xena)
table(is.na(COAD.RPPA.pancan.clean.xena[1,]))
is.na(COAD.RPPA.pancan.clean.xena)
hist(COAD.RPPA.pancan.clean.xena)
COAD.RPPA.pancan.clean.xena.mx <- as.matrix(COAD.RPPA.pancan.clean.xena)
hist(COAD.RPPA.pancan.clean.xena.mx)
boxplot(COAD.RPPA.pancan.clean.xena.mx)

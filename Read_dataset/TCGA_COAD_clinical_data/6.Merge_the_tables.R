##### 5.Merge_the_tables.R
# Merge table TCGA.COAD.Clinical  and TCGA.COAD.Mut
TCGA.COAD.Clinical.output <- dplyr::left_join(TCGA.COAD.Mut, TCGA.COAD.Clinical, by = "patient_barcode")
# Merge table TCGA.COAD.Clinical.output and COAD.SCNA.sample.counts
TCGA.COAD.Clinical.SCAN.output <- dplyr::left_join(TCGA.COAD.Clinical.output, COAD.SCNA.sample.counts, by = "rownames")
table(is.na(TCGA.COAD.Clinical.SCAN.output$SCNA.gene.counts))
saveRDS(TCGA.COAD.Clinical.SCAN.output, file = "2020_2_7_TCGA.COAD.Clinical.Mut.SCNA.output.rds")














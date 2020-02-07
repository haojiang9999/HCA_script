#### 1.1Read_TCGA_dataset.R
COAD_UCSC_Toil_tpm_dataset_merge.syn2623706 <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset_merge.syn2623706.rds")
COAD.pheno <- COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno
COAD.pheno$X_PATIENT
COAD.pheno.merge.syn2623706 <- COAD_UCSC_Toil_tpm_dataset_merge.syn2623706$COAD.pheno.merge.syn2623706

table(COAD.pheno.merge.syn2623706$braf_mut)


clinical_patient_coad$bcr_patient_barcode
table(COAD.pheno$X_PATIENT %in% clinical_patient_coad$bcr_patient_barcode)
colnames(clinical_patient_coad)
table(clinical_patient_coad$microsatellite_instability)
clinical_patient_coad$microsatellite_instability
table(clinical_patient_coad$mismatch_rep_proteins_loss_ihc)
clinical_patient_coad$kras_mutation_found
table(clinical_patient_coad$kras_mutation_found)
clinical_patient_coad[clinical_patient_coad$mismatch_rep_proteins_loss_ihc %in% "CDE_ID:3105496",]
table(clinical_patient_coad$mismatch_rep_proteins_loss_ihc %in% "CDE_ID:3105496")
table(clinical_patient_coad$braf_gene_analysis_result)

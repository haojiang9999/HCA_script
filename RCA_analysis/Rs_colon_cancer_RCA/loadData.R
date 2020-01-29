# load data
#load reference table
referPanel_of_mix <- readRDS("/data8t_4/JH/MyJobs/Colon_SC_Project/Reference_panel/Rs_reference_panel_combine/referPanel_of_mix.rds")
# load expression data set
GSE97693_Tang_TPM_cells <- readRDS("/data8t_4/JH/MyJobs/Colon_SC_Project/GSE97693_Tang_colorectal_cancer/R_RNA_data/GSE97693_Tang_TPM_cells.rds")

BiocManager::install("ComplexHeatmap")
browseVignettes("ComplexHeatmap")
dataSet_for_RCA <- list(referPanel_of_mix = referPanel_of_mix,
                        GSE97693_Tang_TPM_cells = GSE97693_Tang_TPM_cells,
                        GSE97693_Tang_TPM_cells_sample_anno = GSE97693_Tang_TPM_cells_sample_anno)
dataSet_for_RCA
saveRDS(dataSet_for_RCA, file = "dataSet_for_RCA.rds")

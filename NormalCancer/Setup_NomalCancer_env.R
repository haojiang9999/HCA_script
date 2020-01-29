#### setup environment ####
library(clusterProfiler)
### Part1:Cor genes analysis
source("R/Step1_CorValueByGenes_1cell_vs_multiple_reference_cell_types.R")
source("R/Step2_CorValueByGenes_multiple_cells_vs_multiple_reference_cell_types.R")
source("R/Step4_Top_N_genes_in_top_N_cell_types_mx.R")
### Part2:Enrichment analysis
library(devtools)
load_all("clusterProfiler_JH/")
source("R/GeneSymbol2GeneID.R")
source("R/Enrichment_analysis_FUNs_JH_modified.R")
source("R/Merge_enrichment_analysis_resaults.R")
#library()

#install.packages("devtools")

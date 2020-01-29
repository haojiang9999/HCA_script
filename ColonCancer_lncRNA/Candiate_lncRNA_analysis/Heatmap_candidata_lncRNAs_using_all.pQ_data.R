# Using all pQ data
#### Heatmap plot of candidata lncRNA

# construct matrix for plotting
# GeneNames 
geneNmae <- as.character(epi_stem_DE_TvN.lnc.NameFiltered_fc_0.2$X)
### read RCA_all_FPKM_dataset.rds
RCA_all_FPKM_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_all_FPKM_dataset.rds")
NM.all.phenoType <- RCA_all_FPKM_dataset$NM.all.phenoType
tumor.all.phenoType <- RCA_all_FPKM_dataset$tumor.all.phenoType
# expression data
# Find the cell ID in all.data according epithelial.fpkm.exp
NM.stemTA.id <- as.character(NM.Epi.phenoType[samNM.stemTA,]$NM_ID)
tumor.stemTA.id <- as.character(tumor.Epi.phenoType[samTumor.stemTA,]$tumor_ID)
allName.NM <- rownames(NM.all.phenoType)[as.character(NM.all.phenoType$NM_ID) %in% NM.stemTA.id]
allName.tumor <- rownames(tumor.all.phenoType)[as.character(tumor.all.phenoType$tumor_ID) %in% tumor.stemTA.id]

# expression data using all.pQ normalized marix
NM.fpkm.exp.stemTA <- RCA_all.fpkm.exp.pQ[geneNmae, allName.NM]
head(NM.fpkm.exp.stemTA)
# sample names had value in expression data
allName.tumor <- allName.tumor[allName.tumor %in% colnames(RCA_all.fpkm.exp.pQ)]
tumor.fpkm.exp.stemTA <- RCA_all.fpkm.exp.pQ[geneNmae, allName.tumor]
NM.num <- length(allName.NM)
tumor.num <- length(allName.tumor)
data_subset <- as.matrix(cbind(NM.fpkm.exp.stemTA,tumor.fpkm.exp.stemTA))
#### Gene UP or DOWN regulation annotation
table(sign(epi_stem_DE_TvN.lnc.NameFiltered_fc_0.2$DE_fc_tumor_over_normal))
## "Up-tumor" had 26 lncRNAs  "Up-tumor" had 17 lncRNAs
my_gene_row <- data.frame(Gene = as.character(rep(c("Up in tumor", "Down in tumor"), c(26,17))))
rownames(my_gene_row) <- epi_stem_DE_TvN.lnc.NameFiltered_fc_0.2$X
# create heatmap using pheatmap
my_sample_col <- data.frame(sample = rep(c("Normal", "Tumor"), c(NM.num,tumor.num)))
row.names(my_sample_col) <- colnames(data_subset)
data_subset_scaled <- scale(t((data_subset)))
data_subset_scaled <- scale(t(data_subset_scaled))
#data_subset_scaled <- scale(t(log2(data_subset)))
#data_subset_scaled <- scale(t(data_subset_scaled))
library(pheatmap)
pheatmap(data_subset_scaled,annotation_col = my_sample_col,
         annotation_row = my_gene_row,
         scale = c("none"),
         cluster_cols = F, cluster_rows = F)
#### because single sclncRNA data had to see perfect up or down regulated heatmap
# because only part of the cells had pecific lncRNA expressed so 
# I will plot heatmap part by part
### Part 1: Tumor cells
## Up in tumor
Candidate_lncRNA_tumor_up <- pheatmap(data_subset_scaled[1:26,NM.num+1:tumor.num],
                                        cluster_cols = T,cluster_rows = T ,
                                        show_rownames = F, show_colnames = F)
tree_row_tumor_up <- Candidate_lncRNA_tumor_up$tree_row
tree_col_tumor_up <- Candidate_lncRNA_tumor_up$tree_col
## Down in tumor
Candidate_lncRNA_tumor_down <- pheatmap(data_subset_scaled[27:43,NM.num+1:tumor.num],
                                        cluster_cols = tree_col_tumor_up,cluster_rows = T ,
                                        show_rownames = F, show_colnames = F)
tree_row_tumor_down <- Candidate_lncRNA_tumor_down$tree_row
### Part 2: Normal cells
## Down in tumor genes
Candidate_lncRNA_Normal_down <- pheatmap(data_subset_scaled[27:43,1:NM.num],
                                         cluster_cols = T,cluster_rows = tree_row_tumor_down,
                                         show_rownames = T, show_colnames = F)
tree_col_Normal_down <- Candidate_lncRNA_Normal_down$tree_col
## Up in tumor genes
Candidate_lncRNA_Normal_up <- pheatmap(data_subset_scaled[1:26,1:NM.num],
                                         cluster_cols = tree_col_Normal_down,cluster_rows = tree_row_tumor_up,
                                         show_rownames = T, show_colnames = F)


### Save the heatmap
save_pheatmap_png <- function(x, filename, width=2400, height=1000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(Candidate_lncRNA_tumor_up, filename = "Candidate_lncRNA_tumor_upInTu.png")
save_pheatmap_png(Candidate_lncRNA_tumor_down, filename = "Candidate_lncRNA_tumor_downInTu.png")
save_pheatmap_png(Candidate_lncRNA_Normal_up, filename = "Candidate_lncRNA_NM_upInTu.png")
save_pheatmap_png(Candidate_lncRNA_Normal_down, filename = "Candidate_lncRNA_NM_downInTu.png")

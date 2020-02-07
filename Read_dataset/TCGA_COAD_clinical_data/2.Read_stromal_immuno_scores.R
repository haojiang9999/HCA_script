### 2.Read_stromal_immuno_scores
## The Paper:Inferring tumour purity and stromal and immune cell admixture from expression data

ESTIMATE_scores.RNAseq <- read.csv("/data8t_4/JH/MyJobs/Read_dataset/TCGA_COAD_clinical_data/A_list_of_stromal_immune_and_ESTIMATE_scores_in_TCGA_datasets/Data_2_A_list_of_stromal_immune_and_ESTIMATE_scores_in_TCGA_datasets_RNAseq.csv")
ESTIMATE_scores.RNAseqV2 <- read.csv("/data8t_4/JH/MyJobs/Read_dataset/TCGA_COAD_clinical_data/A_list_of_stromal_immune_and_ESTIMATE_scores_in_TCGA_datasets/Data_2_A_list_of_stromal_immune_and_ESTIMATE_scores_in_TCGA_datasets_RNAseqV2.csv")
table(ESTIMATE_scores.RNAseq$ID %in% ESTIMATE_scores.RNAseqV2$ID)

table(ESTIMATE_scores.RNAseq$Platform)
COAD.ESTIMATE.scores.RNAseq <- ESTIMATE_scores.RNAseq[as.character(ESTIMATE_scores.RNAseq$Platform) == "colorectal adenocarcinoma",]
COAD.ESTIMATE.scores.RNAseqV2 <- ESTIMATE_scores.RNAseq[as.character(ESTIMATE_scores.RNAseqV2$Platform) == "colorectal adenocarcinoma",]
table(COAD.ESTIMATE.scores.RNAseqV2$ID %in% COAD.ESTIMATE.scores.RNAseq$ID)
rownames(COAD.pheno.merge.syn2623706)
table(rownames(COAD.pheno.merge.syn2623706) %in% COAD.ESTIMATE.scores.RNAseqV2$ID)
grep("TCGA.A6.2672.01",rownames(COAD.pheno.merge.syn2623706))
grep("TCGA.A6.2672.01",COAD.ESTIMATE.scores.RNAseqV2$ID)
COAD.pheno.merge.syn2623706[rownames(COAD.pheno.merge.syn2623706) %in% ESTIMATE_scores.RNAseqV2$ID,]
ESTIMATE_scores.RNAseqV2[ ESTIMATE_scores.RNAseqV2$ID %in% rownames(COAD.pheno.merge.syn2623706),]


TCGA.A3.3380.01 %in% rownames(COAD.pheno.merge.syn2623706)
cbind(rownames(COAD.pheno.merge.syn2623706),as.character(ESTIMATE_scores.RNAseq$ID))

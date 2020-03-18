#### 3.Survival_data_analysis.R
source("/data8t_4/JH/MyJobs/1_R_script/TCGA_plot/TCGAClusterSurv.R")
# Paper:An Integrated TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR) to drive high quality survival outcome analytics
## 1.Read data
COAD_Survival_SupplementalTabl_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/TCGA_Pan_Cancer/TCGA_Survival_data/COAD_Survival_SupplementalTabl_dataset.rds")
#COAD_Survival_SupplementalTabl_dataset$COAD.Survival_SupplementalTable.metadata
COAD.Survival_SupplementalTable.xena <- COAD_Survival_SupplementalTabl_dataset$COAD.Survival_SupplementalTable.xena
## 2.Merge table
MergeTable.survival <- dplyr::left_join(Cluster.df, COAD.Survival_SupplementalTable.xena, by = "rownames")
## 3.Survival plotting
require(survminer)
library("survival")
require("survival")
fitOS <- survfit(Surv(OS.time, OS) ~ dynamicColors, data = MergeTable.survival)
fitDSS <- survfit(Surv(DSS.time, DSS) ~ dynamicColors, data = MergeTable.survival)
fitDFI <- survfit(Surv(DFI.time, DFI) ~ dynamicColors, data = MergeTable.survival)
fitPFI <- survfit(Surv(PFI.time, PFI) ~ dynamicColors, data = MergeTable.survival)
# Drawing curves
color <- as.character(unique(MergeTable.survival$dynamicColors))
ggsurvplot(fitOS,pval = T,pval.method = T,title = "Over all survival(OS)",  palette =sort(color))

title <- c("Over all survival(OS)", "Disease-specific survival (DSS)",
           "disease-free interval (DFI)","Progression-free interval (PFI)")
ggsurvplot_list(list(fitOS,fitDSS,fitDFI,fitPFI), MergeTable.survival, title = title, risk.table=F,
                pval = T,pval.method = T,
                palette = sort(color))

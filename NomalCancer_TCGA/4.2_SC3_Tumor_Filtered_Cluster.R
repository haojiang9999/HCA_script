#### Tumor cells Cluster using SC3 Filtered
### Step1.Feature selection
# CV filter
Cor.tumor <- Min_max_norm.1500[,TumorID, drop = F]

source("/data8t_4/JH/MyJobs/1_R_script/FUN_TopCV.R")
Cor.tumor.cv10 <- TopCV(Cor.tumor, TopN = 10)
Cor.tumor.cv20 <- TopCV(Cor.tumor, TopN = 20)
Cor.tumor.cv30 <- TopCV(Cor.tumor, TopN = 30)
Cor.tumor.cv40 <- TopCV(Cor.tumor, TopN = 40)
### Step2.Cluster cells (In my script the cells with 0 variance was removed)
SC3.Tumor.cv10 <- JH_SC3_cluster(Cor.tumor.cv10,Pheno.merged.tumor,ks=2:4)
SC3.Tumor.cv20 <- JH_SC3_cluster(Cor.tumor.cv20,Pheno.merged.tumor,ks=2:4)
SC3.Tumor.cv30 <- JH_SC3_cluster(Cor.tumor.cv30,Pheno.merged.tumor,ks=2:4)
SC3.Tumor.cv40 <- JH_SC3_cluster(Cor.tumor.cv40,Pheno.merged.tumor,ks=2:4)
### Step3.Cluster results extraction
## Select a cluster results
hcTumor.cv10 <- SC3.Tumor.cv10$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Cor.tumor.cv10,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hcTumor.cv10)
heatmap.JH(Cor.tumor.cv10,show_colnames = F,
           annotation_col = Pheno.merged.tumor)

hcTumor.cv20 <- SC3.Tumor.cv20$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Cor.tumor.cv20,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hcTumor.cv20)
heatmap.JH(Cor.tumor.cv20,show_colnames = F,
           annotation_col = Pheno.merged.tumor)


hcTumor.cv30 <- SC3.Tumor.cv30$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Cor.tumor.cv30,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hcTumor.cv30)
heatmap.JH(Cor.tumor.cv30,show_colnames = F,
           annotation_col = Pheno.merged.tumor)


hcTumor.cv40 <- SC3.Tumor.cv40$`2`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Cor.tumor.cv40,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hcTumor.cv40)
heatmap.JH(Cor.tumor.cv40,show_colnames = F,
           annotation_col = Pheno.merged.tumor)





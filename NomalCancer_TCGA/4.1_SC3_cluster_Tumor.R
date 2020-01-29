### only tumor samples
TumorID <- rownames(COAD.pheno[COAD.pheno$sampleTypes == "Tumor",])
TumorID
#### Cluster using SC3
### Step1.Load script
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/JH_SC3_cluster.R")
### Step2.Convert NA in matrix to 0
#Cor.tumor <- Min_max_norm.1500[,TumorID, drop = F]
Cor.tumor <- Cor.Res.CV1500$Cor.merged
Cor.tumor<- Cor.tumor[,TumorID]
#Cor.tumor[is.na(Cor.Norm)] <-0
#colnames(Cor.tumor)[colSums(is.na(Cor.tumor)) > 0]
colnames(COAD.pheno)
Pheno.merged.tumor <- COAD.pheno[TumorID,c("sampleTypes","ajcc_pathologic_tumor_stage","histological_type")]

### Step3.Cluster cells (In my script the cells with 0 variance was removed)
SC3.Tumor <- JH_SC3_cluster(Cor.tumor,Pheno.merged.tumor,ks=2:4)
## Select a cluster results
hc1500Tumor <- SC3.Tumor$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
heatmap.JH(Cor.tumor,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hc1500Tumor)


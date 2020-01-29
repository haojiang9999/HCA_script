## 1.Read log10(x+1) transformed scReference 

scReference.log10.CV <- readRDS("/data8t_4/JH/MyJobs/Normal_cell_reference/Step2_Merge_and_Filter_the_scReference/2020_1_19_scReferenceV4.log10.CV.rds")
scReference.list.log10 <- scReference.log10.CV$scReference.list.log10
scReference.list.log10.CV1000 <- scReference.log10.CV$scReference.list.log10.CV1000
scReference.list.log10.CV1500 <- scReference.log10.CV$scReference.list.log10.CV1500
scReference.list.log10.CV2000 <- scReference.log10.CV$scReference.list.log10.CV2000
scReference.list.log10.CV4000 <- scReference.log10.CV$scReference.list.log10.CV4000
scReference.list.log10.CV8000 <- scReference.log10.CV$scReference.list.log10.CV8000


#### 2.Read scRNA colon cancer
### 1).Tang colon cancer
Tang_scColon_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/Tang_scColon_cancer_dataset/Tang_scColon_dataset.rds")
Tang.colon.cancer.FPKM.500 <- Tang_scColon_dataset$FPKM.500.exp
Tang.colon.cancer.TPM.688 <- Tang_scColon_dataset$TPM.688.exp
Tang.colon.cancer.FPKM.500.pheno <- Tang_scColon_dataset$FPKM.500.pheno
Tang.colon.cancer.TPM.688.pheno <- Tang_scColon_dataset$TPM.688.pheno
### 2).RCA colon cancer
RCA_epithelial_COUNT_normal_uniGene_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/GSE81861_RCA_colorectal_tumors/RCA_epithelial_COUNT_normal_uniGene_dataset.rds")
NM.epithelial.COUNT.normal.uniGene <- RCA_epithelial_COUNT_normal_uniGene_dataset$NM.epithelial.COUNT.normal.uniGene
Tumor.epithelial.COUNT.normal.uniGene <- RCA_epithelial_COUNT_normal_uniGene_dataset$Tumor.epithelial.COUNT.normal.uniGene
NM.Epi.phenoType <- RCA_epithelial_COUNT_normal_uniGene_dataset$NM.Epi.phenoType
Tumor.Epi.phenoType <- RCA_epithelial_COUNT_normal_uniGene_dataset$tumor.Epi.phenoType
#### 3.Data transformation 
#### Data transformation ####
#### log10 (x + 1)
log10.Tang.colon.cancer.FPKM.500 <- log10(Tang.colon.cancer.FPKM.500 + 1)
log10.Tang.colon.cancer.TPM.688 <- log10(Tang.colon.cancer.TPM.688 + 1)
log10.NM.epithelial.COUNT.normal.uniGene <- log10(NM.epithelial.COUNT.normal.uniGene + 1)
log10.Tumor.epithelial.COUNT.normal.uniGene <- log10(Tumor.epithelial.COUNT.normal.uniGene + 1)
## build a list
Log10.expList <- list(log10.Tang.colon.cancer.FPKM.500 = log10.Tang.colon.cancer.FPKM.500,
                      log10.Tang.colon.cancer.TPM.688 = log10.Tang.colon.cancer.TPM.688,
                      log10.NM.epithelial.COUNT.normal.uniGene = log10.NM.epithelial.COUNT.normal.uniGene,
                      log10.Tumor.epithelial.COUNT.normal.uniGene = log10.Tumor.epithelial.COUNT.normal.uniGene)
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
##### 4.Distance calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
Cor.Res <- refCorMerge(Log10.expList, scReference.list.log10)
Cor.Res.CV1000 <- refCorMerge(Log10.expList, scReference.list.log10.CV1000)
Cor.Res.CV1500 <- refCorMerge(Log10.expList, scReference.list.log10.CV1500)
Cor.Res.CV2000 <- refCorMerge(Log10.expList, scReference.list.log10.CV2000)
Cor.Res.CV4000 <- refCorMerge(Log10.expList, scReference.list.log10.CV4000)
Cor.Res.CV8000 <- refCorMerge(Log10.expList, scReference.list.log10.CV8000)
Cor.Res.list <- list(Cor.Res = Cor.Res,
                     Cor.Res.CV1000 = Cor.Res.CV1000,
                    Cor.Res.CV1500 = Cor.Res.CV1500,
                    Cor.Res.CV2000 = Cor.Res.CV2000,
                    Cor.Res.CV4000 = Cor.Res.CV4000,
                    Cor.Res.CV8000 = Cor.Res.CV8000)
saveRDS(Cor.Res.list, file = "scColon_log10_screferenceV4_Cor.Res.list.rds")

x <- readRDS("scColon_Pheno.merged.rds")
## cv8000 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/JH_SC3_cluster.R")
Cor.tumor <- Cor.Res.CV8000$Cor.merged
Pheno.merged.tumor <- Pheno.merged
Pheno.merged$
Pheno.merged[,c("Cell_info","CellType")]
### Step3.Cluster cells (In my script the cells with 0 variance was removed)
SC3.Tumor <- JH_SC3_cluster(Cor.tumor,Pheno.merged.tumor,ks=2:4)
## Select a cluster results
hcTumor <- SC3.Tumor$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
options(repr.plot.width=20, repr.plot.height=10)
heatmap.JH(Cor.tumor,show_colnames = F,
           annotation_col = Pheno.merged.tumor, cluster_cols = hcTumor)


## cv2000 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/JH_SC3_cluster.R")
Cor.tumor <- Cor.Res.CV2000$Cor.merged
table(is.na(Cor.tumor))
Cor.tumor[is.na(Cor.tumor)]
colnames(Cor.tumor)[colSums(is.na(Cor.tumor)) > 0]
cell.cor.NA<-Cor.tumor[,c("GSM3241981_scTrioSeq2Rna_CRC09_PT4_546","GSM3443129_scTrioSeq2Rna_CRC09_NC_211")]
View(cell.cor.NA)

Pheno.merged.tumor <- Pheno.merged
### Step3.Cluster cells (In my script the cells with 0 variance was removed)
SC3.Tumor <- JH_SC3_cluster(Cor.tumor,Pheno.merged.tumor,ks=2:4)
## Select a cluster results
hcTumor <- SC3.Tumor$`3`$hc
## Check cluster result
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/heatmap.JH.R")
options(repr.plot.width=20, repr.plot.height=10)
heatmap.JH(Cor.tumor,show_colnames = F,
           annotation_col = Pheno.merged[,c("Cell_info","CellType")], cluster_cols = hcTumor)
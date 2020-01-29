#### Tutorial_CONICSmat_Oligodendroglioma
library(CONICSmat)
devtools::load_all("/data8t_4/JH/MyJobs/test/CONICS/CONICS-master/CONICSmat/")
#### expression table loading
datasetPath <- "/data8t_4/JH/MyJobs/NormalCancer/test_data/dataSet_for_RCA.rds"
Tang_TPM_688 <- readRDS(datasetPath)
Tang_TPM_688.exp <- Tang_TPM_688$GSE97693_Tang_TPM_cells
Tang_TPM_688.anno <- Tang_TPM_688$GSE97693_Tang_TPM_cells_sample_anno
summary(as.numeric(Tang_TPM_688.exp[100,]))
hist(as.numeric(Tang_TPM_688.exp[100,]))
suva_expr = as.matrix(Tang_TPM_688.exp)
## log2(TPM/10+1) expression data 
suva_expr <- log2(suva_expr+0.001)
hist(as.numeric(suva_expr[120,]))
## chang expression <0 to 0
suva_expr [suva_expr <0 ]=0
suva_expr[1:5,1:5]
hist(suva_expr[100,]) # expression data was transformed
## 
dim(suva_expr)
## download matrix containing the coordinates of the human chromosome arms
#download.file(url = "https://raw.githubusercontent.com/diazlab/CONICS/master/chromosome_arm_positions_grch38.txt",
#              destfile = "chromosome_arm_positions_grch38.txt")
regions=read.table("/data8t_4/JH/MyJobs/test/CONICS/chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
head(regions,n=5)
#### Sample dataset: Processing
### Step1 Obtain the chromosomal positions of genes in the expression matrix
gene_pos=getGenePositions(rownames(suva_expr))
### Step2 Subsequently, we can filter uniformative genes. 
suva_expr=filterMatrix(suva_expr,gene_pos[,"hgnc_symbol"],minCells=5)
### Step3 Calculate a normalization factor for each cell.
normFactor=calcNormFactors(suva_expr)
ref=rowMeans(mat[,colnames(mat)[normal]])
normFactor=calcNormFactors(suva_expr)
### Step4 To determine if the average gene 
### expression any of the regions show a bimodal distribution across cells
l=plotAll(suva_expr,normFactor,regions,gene_pos,"SUVA_CNVs")


#### Step5 Plot using my code ####
## Input data
# expression table

### clusterRes data ####
# read AnnoRes for 688
AnnoRes.688.clusterResault <- readRDS(file = "/data8t_4/JH/MyJobs/NormalCancer_688_analysis/AnnoRes.688.clusterResault.rds")
samples <- AnnoRes.688.clusterResault[[1]]$labels
samples.group <- AnnoRes.688.clusterResault[[2]]$groupLabel
samples.anno.df <- data.frame(samples = samples, samples.group = samples.group)
rownames(samples.anno.df) <- samples
##### order cells by cluster #####
samples.anno.df.ordered <- samples.anno.df[order(samples.anno.df$samples.group),]
samples.ordered <- samples.anno.df.ordered$samples
# cluster color table for bar plot
col.df <- samples.anno.df[,2,drop = F]
col.df$samples.group <- sapply(col.df, as.character)
#Pre-set
pmat <- l
expmat <- suva_expr
zscoreThreshold=4
patients=NULL
celltypes=NULL
clusters=2
t=zscoreThreshold
pmat=scale(pmat)
# scale the value
if (max(pmat)>t){
pmat[which(pmat>t)]=t
pmat[which(pmat<(-t))]=(-t)
}else {
mx=min(max(pmat),abs(min(pmat)))
sc=t/mx
pmat=pmat*sc
pmat[which(pmat>t)]=t
pmat[which(pmat<(-t))]=(-t)
}
## correct the colnames
rownames(pmat) <-gsub("_scTrioSeq2Rna_scTrioSeq2Rna_","_scTrioSeq2Rna_",rownames(pmat))
# plot
m <- t(pmat)[,samples.ordered]
pheatmap::pheatmap(m,cluster_rows=F, cluster_cols = F,
                   #cutree_cols = clusters,
                   #show_colnames = F,
                   annotation_col = col.df,
                   col=squash::bluered(100),gaps_col=50,show_colnames = F
                   )




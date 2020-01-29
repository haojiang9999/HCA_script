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





### Step5 Plot begin
##  visualize a heatmap of the posterior probabilities of cells for component2 of each region. 
hi=plotHistogram(l,suva_expr,clusters=6,zscoreThreshold=4,
                 patients=Tang_TPM_688.anno[,"pateintID"],
                 celltypes=Tang_TPM_688.anno[,"cellSitesPatient"] )
table(hi)
## CNV clusters are related to transcriptional signatures
vg=detectVarGenes(suva_expr,500)
#ts=calculateTsne(suva_expr,vg)
#plotTsneGene(ts,suva_expr,c("MBP","CSF1R","ALDOC","OLIG1"))
#plotTsneProbabilities(ts,suva_expr,l[,"1p"],"Tsne Chr1p")
## To obtain a final assignmant as malignant or non-malignant cells,
#lrbic=read.table("SUVA_CNVs_BIC_LR.txt",sep="\t",header=T,row.names=1,check.names=F)
#colnames(lrbic)
#candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & lrbic[,"LRT adj. p-val"]<0.01)]
#hi=plotHistogram(l[,candRegions],suva_expr,clusters=4,zscoreThreshold=4,patients)
## now assign a label as malignant or non-malignant to each cell 
table(Tang_TPM_688.anno[,"cellSites"])
normal= which(Tang_TPM_688.anno[,"cellSites"] == "NC")
tumor=which(Tang_TPM_688.anno[,"cellSites"] != "NC")
## plot the posterior probabilities again, but with statistics for normal and tumor cells:
#redu=plotAll(suva_expr,normFactor,regions[candRegions,],gene_pos,"SUVA_CNVs_with_info.pdf",normal=normal,tumor=tumor)




## visualization of chromosomal alterations in each single cell across the genome for each patient
#par(mfrow=c(1,1))
#plotChromosomeHeatmap(suva_expr,normal = normal, gene_pos = gene_pos, windowsize = 121, chr=T, expThresh=0.2, thresh = 1)
plotChromosomeHeatmap(suva_expr,normal = normal, plotcells = c(normal,tumor) ,gene_pos = gene_pos,chr=T,windowsize = 100, expThresh=0.2, thresh = 1)
plotChromosomeHeatmap()

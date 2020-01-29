#### Tutorial_CONICSmat_Oligodendroglioma
library(CONICSmat)
devtools::load_all("CONICS-master/CONICSmat/")
## log2(TPM/10+1) expression data 
suva_expr = as.matrix(read.table("OG_processed_data_portal.txt",sep="\t",header=T,row.names=1,check.names=F))
suva_expr [which(is.na(suva_expr ))]=0
suva_expr[1:5,1:5]
hist(suva_expr[100,]) # expression data was transformed
## 
dim(suva_expr)
##  information from which patient the cells were derived from
patients=unlist(strsplit(colnames(suva_expr),"_",fixed=TRUE))[seq(1,(3*ncol(suva_expr))-1,3)]
unique(patients)
# To correct the assignment
patients[which(patients=="93")]="MGH93"
patients[which(patients=="97")]="MGH97"

## download matrix containing the coordinates of the human chromosome arms
download.file(url = "https://raw.githubusercontent.com/diazlab/CONICS/master/chromosome_arm_positions_grch38.txt",
              destfile = "chromosome_arm_positions_grch38.txt")
regions=read.table("chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
head(regions,n=5)
#### Sample dataset: Processing
### Step1 Obtain the chromosomal positions of genes in the expression matrix
gene_pos=getGenePositions(rownames(suva_expr))
### Step2 Subsequently, we can filter uniformative genes. 
suva_expr=filterMatrix(suva_expr,gene_pos[,"hgnc_symbol"],minCells=5)
### Step3 Calculate a normalization factor for each cell.
normFactor=calcNormFactors(suva_expr)
### Step4 To determine if the average gene 
### expression any of the regions show a bimodal distribution across cells
l=plotAll(suva_expr,normFactor,regions,gene_pos,"SUVA_CNVs")

### Step5 Plot begin
##  visualize a heatmap of the posterior probabilities of cells for component2 of each region. 
hi=plotHistogram(l,suva_expr,clusters=2,zscoreThreshold=4,patients)
## CNV clusters are related to transcriptional signatures
vg=detectVarGenes(suva_expr,500)
ts=calculateTsne(suva_expr,vg)
plotTsneGene(ts,suva_expr,c("MBP","CSF1R","ALDOC","OLIG1"))
plotTsneProbabilities(ts,suva_expr,l[,"1p"],"Tsne Chr1p")
## To obtain a final assignmant as malignant or non-malignant cells,
lrbic=read.table("SUVA_CNVs_BIC_LR.txt",sep="\t",header=T,row.names=1,check.names=F)
colnames(lrbic)
candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & lrbic[,"LRT adj. p-val"]<0.01)]
hi=plotHistogram(l[,candRegions],suva_expr,clusters=4,zscoreThreshold=4,patients)
## now assign a label as malignant or non-malignant to each cell 
normal= which(hi==1)
tumor=which(hi!=1)
## plot the posterior probabilities again, but with statistics for normal and tumor cells:
redu=plotAll(suva_expr,normFactor,regions[candRegions,],gene_pos,"SUVA_CNVs_with_info.pdf",normal=normal,tumor=tumor)

bin_mat=binarizeMatrix(redu,normal,tumor,0.8)
plotBinaryMat(bin_mat,patients,normal,tumor,patient="MGH97")
detectBreakPoints (suva_expr,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=1,patients=patients,patient="MGH36",breakpoints=regions)
detectBreakPoints (suva_expr,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=4,patients=patients,patient="MGH36",breakpoints=regions)
detectBreakPoints (suva_expr,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=7,patients=patients,patient="MGH36",breakpoints=regions)
detectBreakPoints (suva_expr,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=8,patients=patients,patient="MGH36",breakpoints=regions)
detectBreakPoints (suva_expr,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=10,patients=patients,patient="MGH36",breakpoints=regions)
detectBreakPoints (suva_expr,normal,tumor,windowsize=101,gene_pos=gene_pos,chr=12,patients=patients,patient="MGH36",breakpoints=regions)
## generate a pdf file containing plots for all regions.
plotAllChromosomes (mat = suva_expr,normal = normal,tumor = tumor,windowsize = 101,gene_pos = gene_pos,fname = "MGH36",patients = patients,patient = "MGH36",breakpoints = regions)

## visualization of chromosomal alterations in each single cell across the genome for each patient
par(mfrow=c(1,1))
plotChromosomeHeatmap(suva_expr,normal = normal, plotcells = which(patients=="MGH36"), gene_pos = gene_pos, windowsize = 121, chr=T, expThresh=0.2, thresh = 1)
plotChromosomeHeatmap(suva_expr,normal = normal, plotcells = which(patients=="MGH97") ,gene_pos = gene_pos,chr=T,windowsize = 121, expThresh=0.2, thresh = 1)

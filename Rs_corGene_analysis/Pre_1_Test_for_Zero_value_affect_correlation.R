## load datasets
dataSet_for_RCA <- readRDS("/data8t_4/JH/MyJobs/RCA_analysis/Rs_colon_cancer_RCA/dataSet_for_RCA.rds")
referPanel_of_mix <- dataSet_for_RCA$referPanel_of_mix
GSE97693_Tang_TPM_cells <- dataSet_for_RCA$GSE97693_Tang_TPM_cells
#### Test for zero value affection to two vectors pearson correlation
## Input data two single cells
x <- GSE97693_Tang_TPM_cells[,1]
y <- GSE97693_Tang_TPM_cells[,2]
nonZer0Index <- !x+y == 0 # find zero in both vectors
table(!x+y == 0)          # how many zeros (False)
cov(x,y)
cov(x[nonZer0Index], y[nonZer0Index])
cor(x,y)
cor(x[nonZer0Index], y[nonZer0Index])
### Resaults showed coveriance change a lot, correlation change little

## Input data were one single cells one reference cell type
x <- GSE97693_Tang_TPM_cells[,1,drop = F]
z <- referPanel_of_mix [,1,drop = F]
# find common genes
geneList <- intersect(rownames(x),rownames(z))
i=1
x.com <- x[geneList,,drop = F]
z.com <- z[geneList,,drop = F]
nonZer0Index <- !x.com+z.com == 0
table(nonZer0Index)          # how many zeros (False)
# coveriance and correlation compare
cov(x.com,z.com)
cov(x.com[nonZer0Index], z.com[nonZer0Index])
cor(x.com,z.com)
cor(x.com[nonZer0Index], z.com[nonZer0Index])
##### another problem :
# correlation between reference panel of genes only in single cells and expression > 0  
singleGenesList <- rownames(x.com)[x.com > 0]
cov(x.com,z.com)
cov(x.com[nonZer0Index], z.com[nonZer0Index])
cov(x.com[singleGenesList,], z.com[singleGenesList,])
cor(x.com,z.com)
cor(x.com[nonZer0Index], z.com[nonZer0Index])
cor(x.com[singleGenesList,], z.com[singleGenesList,])

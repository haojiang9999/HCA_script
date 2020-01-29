### find genes have positive coveriance test
x <- GSE97693_Tang_TPM_cells[,1]
y <- GSE97693_Tang_TPM_cells[,2]
## compare 23457 vs 6317 genes 17140 genes were zero in both sample
nonZer0Index <- !x+y == 0 # find zero in both vectors
table(!x+y == 0)          # how many zeros
cov(x,y)
cov(x[nonZer0Index], y[nonZer0Index])
cor(x,y)
cor(x[nonZer0Index], y[nonZer0Index])
## mannully calculate coveriance
# test vector multiple
#x2 <- (x-mean(x))[1:10]
#y2 <- (y-mean(y))[1:10]
#x2 * y2
#34.24133 * -42.63108
sum((x-mean(x))*(y-mean(y)))/length(x)
sum((x[nonZer0Index]-mean(x[nonZer0Index]))*(y[nonZer0Index]-mean(y[nonZer0Index])))/length(nonZer0Index)
## mannully calculate correlation
(sum((x-mean(x))*(y-mean(y)))/length(x))/(sd(x)*sd(y))
(sum((x[nonZer0Index]-mean(x[nonZer0Index]))*(y[nonZer0Index]-mean(y[nonZer0Index])))/length(x[nonZer0Index]))/(sd(x[nonZer0Index])*sd(y[nonZer0Index]))
#### calculate how many genes has positive coveriance
# all genes
covGenes <- (x-mean(x))*(y-mean(y))
which(covGenes < 0)
# expressed genes
covGenes <- (x[nonZer0Index]-mean(x[nonZer0Index]))*(y[nonZer0Index]-mean(y[nonZer0Index]))
covGenes <- covGenes/length(covGenes)

which(covGenes > 0)  # 4839 positive genes
which(covGenes == 0) # 0
which(covGenes < 0)  # 1478 negtive genes
###
#
corValueGenes <- ((x-mean(x))*(y-mean(y))/length(x))/(sd(x)*sd(y))
sort(corValueGenes, decreasing = T)

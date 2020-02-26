#### 2.Significant_Pathway_analysis.R 
### 1.Exame distribution of ssGSEA score
hist(COAD.PanCan33.ssGSEA.xena$TCGA.A6.5657.01)
hist(as.numeric(COAD.PanCan33.ssGSEA.xena[1,]))
COAD.ssGSEA.mx <- as.matrix(COAD.PanCan33.ssGSEA.xena)
plot(density(as.vector(COAD.ssGSEA.mx)), 
     main="COAD.PanCan33.ssGSEA.xena",xlab="ssGSEA score", lwd=2, las=1, xaxt="n", 
     xlim=c(-0.75, 0.75), cex.axis=0.8)
### 2.Using limma to find significant pathway in each groups













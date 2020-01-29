#### sing cell expression analysis between normal and tumor
table(NM.Epi.phenoType$Epi_cellTypes)
table(tumor.Epi.phenoType$Epi_cellTypes)
## stem cell compare
scNM <- rownames(NM.Epi.phenoType)[NM.Epi.phenoType$Epi_cellTypes == "stemTA"]
scTumor <- rownames(tumor.Epi.phenoType)[tumor.Epi.phenoType$Epi_cellTypes == "stemTA"]
#### 1.Data prepare
geneName = "ENSG00000172965"
rowName <- rownames(GeneAnno)[as.character(GeneAnno$GeneID2) == geneName]
## Expression data for normal and cancer
## Using data from pQ normalization
normal.exp<- as.matrix(RCA_epithelial.fpkm.exp.stemTA.pQ[rowName, scNM])
## Samples that in the Normalize matrix
scTumor <- scTumor[scTumor %in% colnames(RCA_epithelial.fpkm.exp.stemTA.pQ)]
tumor.exp <- as.matrix(RCA_epithelial.fpkm.exp.stemTA.pQ[rowName, scTumor])
## fold Change fc
fc <- mean(log2(tumor.exp))/mean(log2(normal.exp))
#fc <- mean(tumor.exp)/mean(normal.exp)
log2(fc)
# Number of tumor and normal samples
NumSam.tumor <- length(tumor.exp)
NumSam.normal <- length(normal.exp)
#### 2.Build ggplot dataframe
df <- data.frame(expression = c(tumor.exp,normal.exp),
                 TissueTypes = c(rep("Tumor",length(tumor.exp)),rep("Normal",length(normal.exp))))
# log2(tpm + 0.01) transform the expression data
df$expressionLog2 <- log2(df$expression)
p_value <- t.test(df$expressionLog2 ~ df$TissueTypes )
#### 3.Arrange ggplot and print the output
## ggplot 
library(ggplot2)
p <- ggplot(df, aes(TissueTypes,expressionLog2))
res <- p + geom_boxplot(fill = c("dodgerblue3","firebrick3")) +  labs(x= paste0("Normal.N = ",NumSam.normal,
                                                                                " Tumor.N = ",NumSam.tumor),
                                                                      y = "Log2(pQ nomalized FPKM) expression",
                                                                      title =paste0(geneName," in ", "RCA_stemTA"),
                                                                      subtitle = paste0("p-value = ", p_value$p.value,
                                                                                        " Fold change = ", fc))
print(res)
ggsave(paste0(geneName,"_in_","RCA_stemTA_scRNA_pQ.pdf"), res,
       width = 6,height = 5)



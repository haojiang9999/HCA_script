#### separete TCGA data into normal and cancer types
TCGA.pheno <- TCGA.Toil.tpm.dataset[[paste0(cancerType,".pheno")]]
## sample type ###
sampleName <- rownames(TCGA.RSEM.gene.tpm)
samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
table(samplType)
## "01" was primary tumor, "11' was normal tissue
tumorIndex <- samplType == "01"
normalIndex <- samplType == "11"
## Expression comparison between normal and cancer
tumor.exp <- TCGA.RSEM.gene.tpm[geneName, tumorIndex]
normal.exp<- TCGA.RSEM.gene.tpm[geneName, normalIndex]
# Number of tumor and normal samples
NumSam.tumor <- length(tumor.exp)
NumSam.normal <- length(normal.exp)
## Build ggplot dataframe
df <- data.frame(expressionLog2_0.01 = c(tumor.exp,normal.exp),
           TissueTypes = c(rep("Tumor",length(tumor.exp)),rep("Normal",length(normal.exp))))
df$expressionLog2_0.01 <- log2(df$expressionLog2_0.01+0.01)
p_value <- t.test(df$expressionLog2_0.01 ~ df$TissueTypes )
p_value$p.value
## ggplot 
library(ggplot2)
p <- ggplot(df, aes(TissueTypes,expressionLog2_0.01))
p + geom_boxplot(fill = c("dodgerblue3","firebrick3")) +  labs(x= paste0("Normal.N = ",NumSam.normal,
                                                                         " Tumor.N = ",NumSam.tumor),
                                                               y = "Log2(tpm + 0.01) expression",
                                                 title =paste0(geneName," in ", cancerType),
                                                 subtitle = paste0("p-value = ", p_value$p.value))

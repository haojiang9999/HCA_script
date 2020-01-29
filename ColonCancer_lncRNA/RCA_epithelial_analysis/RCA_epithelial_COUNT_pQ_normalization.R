#### expression data normalization Using pQ COUNT
library(devtools)
load_all("./NODES_0.0.0.9010/NODES/")
epithelial.COUNT.exp <- cbind(NM.epithelial.COUNT.exp, tumor.epithelial.COUNT.exp)
epithelial.COUNT.exp.pQ <- pQ(data = epithelial.COUNT.exp)
# low expression geen was 31.58095
head(epithelial.COUNT.exp.pQ)

#### DE analysis using NODES
# can not find MetaDE on CRAN
load_all("./MetaDE_1.0.5/MetaDE/")
DE.COUNT <- NODES(data = epithelial.COUNT.exp.pQ, group = c(rep("Normal",length(normal.exp)),rep("tumor",length(tumor.exp))))
geneName = "ENSG00000255198"
rowName <- rownames(GeneAnno)[as.character(GeneAnno$GeneID2) == geneName]
DE.COUNT[rowName,]
DE.fpkm[rowName,]

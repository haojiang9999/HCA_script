#### expression data normalization Using pQ
library(devtools)
load_all("./NODES_0.0.0.9010/NODES/")
epithelial.fpkm.exp.stemTA <- cbind(NM.epithelial.fpkm.exp[,scNM], tumor.epithelial.fpkm.exp[,scTumor])
epithelial.fpkm.exp.stemTA.pQ <- pQ(data = epithelial.fpkm.exp.stemTA, frac = 0.8)
saveRDS()
# low expression geen was 31.58095
head(epithelial.fpkm.exp.stemTA.pQ)

#### DE analysis using NODES
# can not find MetaDE on CRAN
load_all("./MetaDE_1.0.5/MetaDE/")
DE.fpkm <- NODES(data = epithelial.fpkm.exp.stemTA.pQ, group = c(rep("Normal",length(scNM)),rep("tumor",length(scTumor))))
geneName = "ENSG00000248527"
rowName <- rownames(GeneAnno)[as.character(GeneAnno$GeneID2) == geneName]
DE.fpkm[rowName,]
#BiocManager::install("MetaDE")

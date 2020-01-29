#### epithelial.fpkm.exp.stemTA expression data normalization Using pQ
library(devtools)
load_all("./NODES_0.0.0.9010/NODES/")
## stem cell compare
scNM <- rownames(NM.Epi.phenoType)[NM.Epi.phenoType$Epi_cellTypes == "stemTA"]
scTumor <- rownames(tumor.Epi.phenoType)[tumor.Epi.phenoType$Epi_cellTypes == "stemTA"]
epithelial.fpkm.exp.stemTA <- cbind(NM.epithelial.fpkm.exp[,scNM], tumor.epithelial.fpkm.exp[,scTumor])
epithelial.fpkm.exp.stemTA.pQ <- pQ(data = epithelial.fpkm.exp.stemTA)
saveRDS(epithelial.fpkm.exp.stemTA.pQ, file = "RCA_epithelial.fpkm.exp.stemTA.pQ.rds")

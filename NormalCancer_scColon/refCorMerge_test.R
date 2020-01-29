##test of refCorMerge

mxList <- list(Tang.colon.cancer.FPKM.500=Tang.colon.cancer.FPKM.500,
               Tang.colon.cancer.TPM.688 = Tang.colon.cancer.TPM.688,
               NM.epithelial.COUNT.normal.uniGene=NM.epithelial.COUNT.normal.uniGene,
               Tumor.epithelial.COUNT.normal.uniGene = Tumor.epithelial.COUNT.normal.uniGene)
ReferenceList <- log2.scReference.list.CV.8000

Res <- refCorMerge(mxList, ReferenceList)

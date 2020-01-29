#### Correlation of scColon cancer with Normal cancer
source("/data8t_4/JH/MyJobs/1_R_script/RCA_seperate_fuction.R")
source("/data8t_4/JH/MyJobs/1_R_script/refCor.R")
#### 1.Correlation of cancer cells with unfiltered scReference.V1
Cor.Tang.colon.cancer.FPKM.500 <- refCor(Tang.colon.cancer.FPKM.500, scReference.V1[["scReference.list"]])
Cor.Tang.colon.cancer.TPM.688 <- refCor(Tang.colon.cancer.TPM.688, scReference.V1[["scReference.list"]])
Cor.RCA.NM.epithelial.COUNT.normal <- refCor(NM.epithelial.COUNT.normal.uniGene, scReference.V1[["scReference.list"]])
Cor.RCA.Tumor.epithelial.COUNT.normal <- refCor(Tumor.epithelial.COUNT.normal.uniGene, scReference.V1[["scReference.list"]])
Cor.Orgnoid.colon.TPM <- refCor(Orgnoid_expression.uniGene.TPM, scReference.V1[["scReference.list"]])

#### 2.Merge to one table to analysis
Cor.merged<- cbind(Cor.Tang.colon.cancer.FPKM.500,
      Cor.Tang.colon.cancer.TPM.688,
      Cor.RCA.NM.epithelial.COUNT.normal,
      Cor.RCA.Tumor.epithelial.COUNT.normal)


#### 3.Build the sample phenoType infor
Pheno.merged <- as.data.frame(colnames(Cor.merged))
colnames(Pheno.merged) <- "Cell_names"
rownames(Pheno.merged) <- Pheno.merged$Cell_names
Pheno.merged$Cell_info <- c(rep("Cor.Tang.colon.cancer.FPKM.500", ncol(Cor.Tang.colon.cancer.FPKM.500)),
                            rep("Cor.Tang.colon.cancer.TPM.688", ncol(Cor.Tang.colon.cancer.TPM.688)),
                            rep("Cor.RCA.NM.epithelial.COUNT.normal", ncol(Cor.RCA.NM.epithelial.COUNT.normal)),
                            rep("Cor.RCA.Tumor.epithelial.COUNT.normal", ncol(Cor.RCA.Tumor.epithelial.COUNT.normal)))










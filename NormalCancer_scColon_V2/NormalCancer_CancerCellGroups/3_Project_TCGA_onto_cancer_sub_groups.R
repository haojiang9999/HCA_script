#### Project TCGA data onto cancer cell groups
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
mxList <- list(COAD_tpm_symbol =  COAD_tpm_symbol) 
test <- refCorMerge(mxList, Cancer.group.list$G1.exp)

Cor.merged <- test$Cor.merged
head(COAD.RSEM.gene.tpm)
View(Cor.merged)

#### Gene filter version
geneList <- intersect(rownames(COAD_tpm_symbol),genes.2000)
COAD_symbol_sub <- COAD_tpm_symbol[geneList,]
mxList <- list(COAD_symbol_sub =  COAD_symbol_sub) 
Cor_G1 <- refCorMerge(mxList, Cancer.group.list$G1.exp)
Cor_G2 <- refCorMerge(mxList, Cancer.group.list$G2.exp)
Cor_G3 <- refCorMerge(mxList, Cancer.group.list$G3.exp)
Cor_G4 <- refCorMerge(mxList, Cancer.group.list$G4.exp)
Cor_G5 <- refCorMerge(mxList, Cancer.group.list$G5.exp)
Cor_G6 <- refCorMerge(mxList, Cancer.group.list$G6.exp)
Cor.merged <- test$Cor.merged
head(COAD.RSEM.gene.tpm)
View(Cor_G1$Cor.merged)
View(Cor_G2$Cor.merged)
x <- na.omit(Cor_G1$Cor.merged)
summary(apply(x , 2, median))
summary(apply(na.omit(Cor_G2$Cor.merged) , 2, median))
summary(apply(na.omit(Cor_G3$Cor.merged) , 2, median))
summary(apply(na.omit(Cor_G4$Cor.merged) , 2, median))
summary(apply(na.omit(Cor_G5$Cor.merged) , 2, median))

test<- rbind(apply(na.omit(Cor_G1$Cor.merged) , 2, median),
             apply(na.omit(Cor_G2$Cor.merged) , 2, median),
      apply(na.omit(Cor_G3$Cor.merged) , 2, median),
      apply(na.omit(Cor_G4$Cor.merged) , 2, median),
      apply(na.omit(Cor_G5$Cor.merged) , 2, median),
      apply(na.omit(Cor_G6$Cor.merged) , 2, median))
test2<- rbind(apply(na.omit(Cor_G1$Cor.merged) , 2, mean),
             apply(na.omit(Cor_G2$Cor.merged) , 2, mean),
             apply(na.omit(Cor_G3$Cor.merged) , 2, mean),
             apply(na.omit(Cor_G4$Cor.merged) , 2, mean),
             apply(na.omit(Cor_G5$Cor.merged) , 2, mean),
             apply(na.omit(Cor_G6$Cor.merged) , 2, mean))

####### Average the duplicated genes in COUNT_normal_dataset ######
GeneID2 <- GeneAnno$GeneID2
library(clusterProfiler)
eg = bitr(GeneID2, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
GeneAnno.sub <- GeneAnno[GeneAnno$GeneID2 %in%  eg$ENSEMBL,]
table(duplicated(GeneAnno.sub$GeneSymbol))
# Symbol of dup or non-dup genes
Dup.genes <- as.character(unique(GeneAnno.sub[duplicated(GeneAnno.sub$GeneSymbol),]$GeneSymbol))
#NonDup.genes <-as.character(unique(GeneAnno.sub[!duplicated(GeneAnno.sub$GeneSymbol),]$GeneSymbol))
NonDup.genes <-as.character(GeneAnno.sub[!GeneAnno.sub$GeneSymbol %in% Dup.genes,]$GeneSymbol)
### Expression dataframe of non dup genes
Index.NonDup<- rownames(GeneAnno.sub[!GeneAnno.sub$GeneSymbol %in% Dup.genes,])## remove dup genes
NM.epithelial.COUNT.normal.uni1 <- NM.epithelial.COUNT.normal[Index.NonDup,]
Tumor.epithelial.COUNT.normal.uni1 <- Tumor.epithelial.COUNT.normal[Index.NonDup,]
# Change row names to gene symbol
rownames(NM.epithelial.COUNT.normal.uni1) <- GeneAnno.sub[Index.NonDup,]$GeneSymbol
rownames(Tumor.epithelial.COUNT.normal.uni1) <- GeneAnno.sub[Index.NonDup,]$GeneSymbol
### Average the geme expression of duplicated genes
NM.epithelial.COUNT.normal.uni2 <- sapply(Dup.genes, FUN = function(symbol){
  #symbol= "SNORA63"
  Index.Dup<- rownames(GeneAnno.sub[GeneAnno.sub$GeneSymbol %in% symbol,])
  NM.uni2 <- NM.epithelial.COUNT.normal[Index.Dup,]
  #Tumor.uni2 <- Tumor.epithelial.COUNT.normal[Index.Dup,]
  NM.uni2 <- Matrix::colMeans(NM.uni2)
  #Tumor.uni2 <- Matrix::colMeans(Tumor.uni2)
})
NM.epithelial.COUNT.normal.uni2 <- t(NM.epithelial.COUNT.normal.uni2)
Tumor.epithelial.COUNT.normal.uni2 <- sapply(Dup.genes, FUN = function(symbol){
  #symbol= "SNORA63"
  Index.Dup<- rownames(GeneAnno.sub[GeneAnno.sub$GeneSymbol %in% symbol,])
  #NM.uni2 <- NM.epithelial.COUNT.normal[Index.Dup,]
  Tumor.uni2 <- Tumor.epithelial.COUNT.normal[Index.Dup,]
  #NM.uni2 <- Matrix::colMeans(NM.uni2)
  Tumor.uni2 <- Matrix::colMeans(Tumor.uni2)
})
Tumor.epithelial.COUNT.normal.uni2 <- t(Tumor.epithelial.COUNT.normal.uni2)
#### Combine two unique table
NM.epithelial.COUNT.normal.uniGene <- rbind(NM.epithelial.COUNT.normal.uni1,NM.epithelial.COUNT.normal.uni2)
Tumor.epithelial.COUNT.normal.uniGene <- rbind(Tumor.epithelial.COUNT.normal.uni1,Tumor.epithelial.COUNT.normal.uni2)
table(duplicated(rownames(NM.epithelial.COUNT.normal.uniGene)))
table(duplicated(rownames(Tumor.epithelial.COUNT.normal.uniGene)))

# Build dataset 
RCA_epithelial_COUNT_normal_uniGene_dataset <- list(NM.epithelial.COUNT.normal.uniGene = NM.epithelial.COUNT.normal.uniGene,
                                                    Tumor.epithelial.COUNT.normal.uniGene = Tumor.epithelial.COUNT.normal.uniGene,
                                                    NM.Epi.phenoType = NM.Epi.phenoType,
                                                    tumor.Epi.phenoType = tumor.Epi.phenoType,
                                                    GeneAnno = GeneAnno)
saveRDS(RCA_epithelial_COUNT_normal_uniGene_dataset, file = "RCA_epithelial_COUNT_normal_uniGene_dataset.rds")

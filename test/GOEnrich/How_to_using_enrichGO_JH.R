library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
GO_DATA_BP <- readRDS("GO_DATA_BP.rds")
GO_DATA_CC <- readRDS("GO_DATA_CC.rds")
GO_DATA_MF <- readRDS("GO_DATA_MF.rds")

ego_CC <- enrichGO_JH(gene          = gene,
               universe      = names(geneList),
               OrgDb         = org.Hs.eg.db,
               ont           = "CC",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05,
               readable      = TRUE)
head(ego_CC)

ego_BP <- enrichGO_JH(gene          = gene,
                      universe      = names(geneList),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
head(ego_BP)
OrgDb = org.Hs.eg.db
test <- enrichGO_JH(gene= gene,ont = "BP")

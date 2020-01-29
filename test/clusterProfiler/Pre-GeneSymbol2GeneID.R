#### pre-geneID convert.r ####
### gene symbol to geneID
library(clusterProfiler)
library(org.Hs.eg.db)
ids1 = bitr(input, fromType="SYMBOL", toType="ENTREZID", 
            OrgDb="org.Hs.eg.db",drop = F)
head(ids1)
leftSymbol <- ids1[is.na(ids1[,2]),][,1]
ids2 = bitr(leftSymbol, fromType="ALIAS", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(ids2)
colnames(ids2)[1] <- "SYMBOL"
ids <- rbind(ids1,ids2)

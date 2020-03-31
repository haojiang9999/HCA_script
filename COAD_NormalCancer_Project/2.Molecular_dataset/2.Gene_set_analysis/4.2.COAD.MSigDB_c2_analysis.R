#### 4.2.COAD.MSigDB_c2_analysis.R
dt.c2 <- decideTests(fit.c2.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.c2.all <- topTable(fit.c2.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.c2.all <- topTable(fit.c2.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.c2.all <- topTable(fit.c2.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.c2.all <- topTable(fit.c2.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.c2.all <- topTable(fit.c2.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.c2.all <- topTable(fit.c2.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.c2.all)
dim(DEgeneSets.brown.c2.all)
dim(DEgeneSets.turquoise.c2.all)
dim(DEgeneSets.BlvBr.c2.all)

head(DEgeneSets.blue.c2.all)


#cbind(rownames(DEgeneSets.blue.c2.all),rownames(DEgeneSets.brown.c2.all),rownames(DEgeneSets.turquoise.c2.all))
## Order by rownames
Index.c2 <- rownames(DEgeneSets.blue.c2.all)
## Ordered q value table log10transformed q value
c2.log10.q.blue <- -log10(DEgeneSets.blue.c2.all[Index.c2,]$adj.P.Val) * sign(DEgeneSets.blue.c2.all[Index.c2,]$logFC)
c2.log10.q.brown <- -log10(DEgeneSets.brown.c2.all[Index.c2,]$adj.P.Val) * sign(DEgeneSets.brown.c2.all[Index.c2,]$logFC)
c2.log10.q.turquoise <- -log10(DEgeneSets.turquoise.c2.all[Index.c2,]$adj.P.Val) * sign(DEgeneSets.turquoise.c2.all[Index.c2,]$logFC)
c2.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.c2.all[Index.c2,]$adj.P.Val) * sign(DEgeneSets.BlvBr.c2.all[Index.c2,]$logFC)
c2.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.c2.all[Index.c2,]$adj.P.Val) * sign(DEgeneSets.BlvTu.c2.all[Index.c2,]$logFC)
c2.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.c2.all[Index.c2,]$adj.P.Val) * sign(DEgeneSets.BrvTu.c2.all[Index.c2,]$logFC)

c2.log10.q.all <- cbind(c2.log10.q.blue,c2.log10.q.brown,c2.log10.q.turquoise,
                        c2.log10.q.BlvBr,c2.log10.q.BlvTu,c2.log10.q.BrvTu)
rownames(c2.log10.q.all) <- Index.c2
head(c2.log10.q.all)
saveRDS(c2.log10.q.all, file = "COAD.tb.c2.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
c2.log10.q.all.m <- melt(c2.log10.q.all[1:10,])
ggplot(c2.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

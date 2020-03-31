#### 4.5.COAD.MSigDB_c4_analysis.R
dt.c4 <- decideTests(fit.c4.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.c4.all <- topTable(fit.c4.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.c4.all <- topTable(fit.c4.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.c4.all <- topTable(fit.c4.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.c4.all <- topTable(fit.c4.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.c4.all <- topTable(fit.c4.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.c4.all <- topTable(fit.c4.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.c4.all)
dim(DEgeneSets.brown.c4.all)
dim(DEgeneSets.turquoise.c4.all)
dim(DEgeneSets.BlvBr.c4.all)

head(DEgeneSets.blue.c4.all)


#cbind(rownames(DEgeneSets.blue.c4.all),rownames(DEgeneSets.brown.c4.all),rownames(DEgeneSets.turquoise.c4.all))
## Order by rownames
Index.c4 <- rownames(DEgeneSets.blue.c4.all)
## Ordered q value table log10transformed q value
c4.log10.q.blue <- -log10(DEgeneSets.blue.c4.all[Index.c4,]$adj.P.Val) * sign(DEgeneSets.blue.c4.all[Index.c4,]$logFC)
c4.log10.q.brown <- -log10(DEgeneSets.brown.c4.all[Index.c4,]$adj.P.Val) * sign(DEgeneSets.brown.c4.all[Index.c4,]$logFC)
c4.log10.q.turquoise <- -log10(DEgeneSets.turquoise.c4.all[Index.c4,]$adj.P.Val) * sign(DEgeneSets.turquoise.c4.all[Index.c4,]$logFC)
c4.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.c4.all[Index.c4,]$adj.P.Val) * sign(DEgeneSets.BlvBr.c4.all[Index.c4,]$logFC)
c4.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.c4.all[Index.c4,]$adj.P.Val) * sign(DEgeneSets.BlvTu.c4.all[Index.c4,]$logFC)
c4.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.c4.all[Index.c4,]$adj.P.Val) * sign(DEgeneSets.BrvTu.c4.all[Index.c4,]$logFC)

c4.log10.q.all <- cbind(c4.log10.q.blue,c4.log10.q.brown,c4.log10.q.turquoise,
                        c4.log10.q.BlvBr,c4.log10.q.BlvTu,c4.log10.q.BrvTu)
rownames(c4.log10.q.all) <- Index.c4
head(c4.log10.q.all)
saveRDS(c4.log10.q.all, file = "COAD.tb.c4.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
c4.log10.q.all.m <- melt(c4.log10.q.all[1:10,])
ggplot(c4.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()


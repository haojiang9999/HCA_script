#### 4.6.COAD.MSigDB_c5_analysis.R
dt.c5 <- decideTests(fit.c5.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.c5.all <- topTable(fit.c5.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.c5.all <- topTable(fit.c5.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.c5.all <- topTable(fit.c5.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.c5.all <- topTable(fit.c5.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.c5.all <- topTable(fit.c5.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.c5.all <- topTable(fit.c5.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.c5.all)
dim(DEgeneSets.brown.c5.all)
dim(DEgeneSets.turquoise.c5.all)
dim(DEgeneSets.BlvBr.c5.all)

head(DEgeneSets.blue.c5.all)


#cbind(rownames(DEgeneSets.blue.c5.all),rownames(DEgeneSets.brown.c5.all),rownames(DEgeneSets.turquoise.c5.all))
## Order by rownames
Index.c5 <- rownames(DEgeneSets.blue.c5.all)
## Ordered q value table log10transformed q value
c5.log10.q.blue <- -log10(DEgeneSets.blue.c5.all[Index.c5,]$adj.P.Val) * sign(DEgeneSets.blue.c5.all[Index.c5,]$logFC)
c5.log10.q.brown <- -log10(DEgeneSets.brown.c5.all[Index.c5,]$adj.P.Val) * sign(DEgeneSets.brown.c5.all[Index.c5,]$logFC)
c5.log10.q.turquoise <- -log10(DEgeneSets.turquoise.c5.all[Index.c5,]$adj.P.Val) * sign(DEgeneSets.turquoise.c5.all[Index.c5,]$logFC)
c5.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.c5.all[Index.c5,]$adj.P.Val) * sign(DEgeneSets.BlvBr.c5.all[Index.c5,]$logFC)
c5.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.c5.all[Index.c5,]$adj.P.Val) * sign(DEgeneSets.BlvTu.c5.all[Index.c5,]$logFC)
c5.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.c5.all[Index.c5,]$adj.P.Val) * sign(DEgeneSets.BrvTu.c5.all[Index.c5,]$logFC)

c5.log10.q.all <- cbind(c5.log10.q.blue,c5.log10.q.brown,c5.log10.q.turquoise,
                        c5.log10.q.BlvBr,c5.log10.q.BlvTu,c5.log10.q.BrvTu)
rownames(c5.log10.q.all) <- Index.c5
head(c5.log10.q.all)
saveRDS(c5.log10.q.all, file = "COAD.tb.c5.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
c5.log10.q.all.m <- melt(c5.log10.q.all[1:10,])
ggplot(c5.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()


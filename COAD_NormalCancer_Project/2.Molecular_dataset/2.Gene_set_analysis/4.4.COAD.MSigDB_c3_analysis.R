#### 4.4.COAD.MSigDB_c3_analysis.R
dt.c3 <- decideTests(fit.c3.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.c3.all <- topTable(fit.c3.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.c3.all <- topTable(fit.c3.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.c3.all <- topTable(fit.c3.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.c3.all <- topTable(fit.c3.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.c3.all <- topTable(fit.c3.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.c3.all <- topTable(fit.c3.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.c3.all)
dim(DEgeneSets.brown.c3.all)
dim(DEgeneSets.turquoise.c3.all)
dim(DEgeneSets.BlvBr.c3.all)

head(DEgeneSets.blue.c3.all)


#cbind(rownames(DEgeneSets.blue.c3.all),rownames(DEgeneSets.brown.c3.all),rownames(DEgeneSets.turquoise.c3.all))
## Order by rownames
Index.c3 <- rownames(DEgeneSets.blue.c3.all)
## Ordered q value table log10transformed q value
c3.log10.q.blue <- -log10(DEgeneSets.blue.c3.all[Index.c3,]$adj.P.Val) * sign(DEgeneSets.blue.c3.all[Index.c3,]$logFC)
c3.log10.q.brown <- -log10(DEgeneSets.brown.c3.all[Index.c3,]$adj.P.Val) * sign(DEgeneSets.brown.c3.all[Index.c3,]$logFC)
c3.log10.q.turquoise <- -log10(DEgeneSets.turquoise.c3.all[Index.c3,]$adj.P.Val) * sign(DEgeneSets.turquoise.c3.all[Index.c3,]$logFC)
c3.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.c3.all[Index.c3,]$adj.P.Val) * sign(DEgeneSets.BlvBr.c3.all[Index.c3,]$logFC)
c3.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.c3.all[Index.c3,]$adj.P.Val) * sign(DEgeneSets.BlvTu.c3.all[Index.c3,]$logFC)
c3.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.c3.all[Index.c3,]$adj.P.Val) * sign(DEgeneSets.BrvTu.c3.all[Index.c3,]$logFC)

c3.log10.q.all <- cbind(c3.log10.q.blue,c3.log10.q.brown,c3.log10.q.turquoise,
                        c3.log10.q.BlvBr,c3.log10.q.BlvTu,c3.log10.q.BrvTu)
rownames(c3.log10.q.all) <- Index.c3
head(c3.log10.q.all)
saveRDS(c3.log10.q.all, file = "COAD.tb.c3.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
c3.log10.q.all.m <- melt(c3.log10.q.all[1:10,])
ggplot(c3.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()


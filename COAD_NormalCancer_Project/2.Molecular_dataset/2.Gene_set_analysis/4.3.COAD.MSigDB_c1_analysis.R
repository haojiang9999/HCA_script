#### 4.3.COAD.MSigDB_c1_analysis.R
dt.c1 <- decideTests(fit.c1.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.c1.all <- topTable(fit.c1.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.c1.all <- topTable(fit.c1.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.c1.all <- topTable(fit.c1.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.c1.all <- topTable(fit.c1.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.c1.all <- topTable(fit.c1.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.c1.all <- topTable(fit.c1.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.c1.all)
dim(DEgeneSets.brown.c1.all)
dim(DEgeneSets.turquoise.c1.all)
dim(DEgeneSets.BlvBr.c1.all)

head(DEgeneSets.blue.c1.all)


#cbind(rownames(DEgeneSets.blue.c1.all),rownames(DEgeneSets.brown.c1.all),rownames(DEgeneSets.turquoise.c1.all))
## Order by rownames
Index.c1 <- rownames(DEgeneSets.blue.c1.all)
## Ordered q value table log10transformed q value
c1.log10.q.blue <- -log10(DEgeneSets.blue.c1.all[Index.c1,]$adj.P.Val) * sign(DEgeneSets.blue.c1.all[Index.c1,]$logFC)
c1.log10.q.brown <- -log10(DEgeneSets.brown.c1.all[Index.c1,]$adj.P.Val) * sign(DEgeneSets.brown.c1.all[Index.c1,]$logFC)
c1.log10.q.turquoise <- -log10(DEgeneSets.turquoise.c1.all[Index.c1,]$adj.P.Val) * sign(DEgeneSets.turquoise.c1.all[Index.c1,]$logFC)
c1.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.c1.all[Index.c1,]$adj.P.Val) * sign(DEgeneSets.BlvBr.c1.all[Index.c1,]$logFC)
c1.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.c1.all[Index.c1,]$adj.P.Val) * sign(DEgeneSets.BlvTu.c1.all[Index.c1,]$logFC)
c1.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.c1.all[Index.c1,]$adj.P.Val) * sign(DEgeneSets.BrvTu.c1.all[Index.c1,]$logFC)

c1.log10.q.all <- cbind(c1.log10.q.blue,c1.log10.q.brown,c1.log10.q.turquoise,
                        c1.log10.q.BlvBr,c1.log10.q.BlvTu,c1.log10.q.BrvTu)
rownames(c1.log10.q.all) <- Index.c1
head(c1.log10.q.all)
saveRDS(c1.log10.q.all, file = "COAD.tb.c1.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
c1.log10.q.all.m <- melt(c1.log10.q.all[1:10,])
ggplot(c1.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()


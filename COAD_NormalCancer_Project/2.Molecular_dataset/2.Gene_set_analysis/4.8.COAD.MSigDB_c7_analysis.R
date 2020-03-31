#### 4.8.COAD.MSigDB_c7_analysis.R
dt.c7 <- decideTests(fit.c7.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.c7.all <- topTable(fit.c7.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.c7.all <- topTable(fit.c7.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.c7.all <- topTable(fit.c7.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.c7.all <- topTable(fit.c7.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.c7.all <- topTable(fit.c7.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.c7.all <- topTable(fit.c7.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.c7.all)
dim(DEgeneSets.brown.c7.all)
dim(DEgeneSets.turquoise.c7.all)
dim(DEgeneSets.BlvBr.c7.all)

head(DEgeneSets.blue.c7.all)


#cbind(rownames(DEgeneSets.blue.c7.all),rownames(DEgeneSets.brown.c7.all),rownames(DEgeneSets.turquoise.c7.all))
## Order by rownames
Index.c7 <- rownames(DEgeneSets.blue.c7.all)
## Ordered q value table log10transformed q value
c7.log10.q.blue <- -log10(DEgeneSets.blue.c7.all[Index.c7,]$adj.P.Val) * sign(DEgeneSets.blue.c7.all[Index.c7,]$logFC)
c7.log10.q.brown <- -log10(DEgeneSets.brown.c7.all[Index.c7,]$adj.P.Val) * sign(DEgeneSets.brown.c7.all[Index.c7,]$logFC)
c7.log10.q.turquoise <- -log10(DEgeneSets.turquoise.c7.all[Index.c7,]$adj.P.Val) * sign(DEgeneSets.turquoise.c7.all[Index.c7,]$logFC)
c7.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.c7.all[Index.c7,]$adj.P.Val) * sign(DEgeneSets.BlvBr.c7.all[Index.c7,]$logFC)
c7.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.c7.all[Index.c7,]$adj.P.Val) * sign(DEgeneSets.BlvTu.c7.all[Index.c7,]$logFC)
c7.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.c7.all[Index.c7,]$adj.P.Val) * sign(DEgeneSets.BrvTu.c7.all[Index.c7,]$logFC)

c7.log10.q.all <- cbind(c7.log10.q.blue,c7.log10.q.brown,c7.log10.q.turquoise,
                        c7.log10.q.BlvBr,c7.log10.q.BlvTu,c7.log10.q.BrvTu)
rownames(c7.log10.q.all) <- Index.c7
head(c7.log10.q.all)
saveRDS(c7.log10.q.all, file = "COAD.tb.c7.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
c7.log10.q.all.m <- melt(c7.log10.q.all[1:10,])
ggplot(c7.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()


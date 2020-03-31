#### 4.1.COAD.MSigDB_h_analysis.R
dt.h <- decideTests(fit.h.2, p.value = adjPvalueCutoff )
#number = Inf
DEgeneSets.blue.h.all <- topTable(fit.h.2, coef="blue - (brown + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.brown.h.all <- topTable(fit.h.2, coef="brown - (blue + turquoise)/2", number=Inf,adjust="BH")
DEgeneSets.turquoise.h.all <- topTable(fit.h.2, coef="turquoise - (brown + blue)/2", number=Inf,adjust="BH")
DEgeneSets.BlvBr.h.all <- topTable(fit.h.2, coef="blue - brown",number=Inf,adjust="BH")
DEgeneSets.BlvTu.h.all <- topTable(fit.h.2, coef="blue - turquoise", number=Inf,adjust="BH")
DEgeneSets.BrvTu.h.all <- topTable(fit.h.2, coef="brown - turquoise", number=Inf,adjust="BH")
dim(DEgeneSets.blue.h.all)
dim(DEgeneSets.brown.h.all)
dim(DEgeneSets.turquoise.h.all)
dim(DEgeneSets.BlvBr.h.all)

head(DEgeneSets.blue.h.all)


#cbind(rownames(DEgeneSets.blue.h.all),rownames(DEgeneSets.brown.h.all),rownames(DEgeneSets.turquoise.h.all))
## Order by rownames
Index.h <- rownames(DEgeneSets.blue.h.all)
## Ordered q value table log10transformed q value
h.log10.q.blue <- -log10(DEgeneSets.blue.h.all[Index.h,]$adj.P.Val) * sign(DEgeneSets.blue.h.all[Index.h,]$logFC)
h.log10.q.brown <- -log10(DEgeneSets.brown.h.all[Index.h,]$adj.P.Val) * sign(DEgeneSets.brown.h.all[Index.h,]$logFC)
h.log10.q.turquoise <- -log10(DEgeneSets.turquoise.h.all[Index.h,]$adj.P.Val) * sign(DEgeneSets.turquoise.h.all[Index.h,]$logFC)
h.log10.q.BlvBr <- -log10(DEgeneSets.BlvBr.h.all[Index.h,]$adj.P.Val) * sign(DEgeneSets.BlvBr.h.all[Index.h,]$logFC)
h.log10.q.BlvTu <- -log10(DEgeneSets.BlvTu.h.all[Index.h,]$adj.P.Val) * sign(DEgeneSets.BlvTu.h.all[Index.h,]$logFC)
h.log10.q.BrvTu <- -log10(DEgeneSets.BrvTu.h.all[Index.h,]$adj.P.Val) * sign(DEgeneSets.BrvTu.h.all[Index.h,]$logFC)

h.log10.q.all <- cbind(h.log10.q.blue,h.log10.q.brown,h.log10.q.turquoise,
                        h.log10.q.BlvBr,h.log10.q.BlvTu,h.log10.q.BrvTu)
rownames(h.log10.q.all) <- Index.h
head(h.log10.q.all)
saveRDS(h.log10.q.all, file = "COAD.tb.h.log10.q.all.rds")

### Heatmap plot example
library(reshape2)
library(ggplot2)
h.log10.q.all.m <- melt(h.log10.q.all)
ggplot(h.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()



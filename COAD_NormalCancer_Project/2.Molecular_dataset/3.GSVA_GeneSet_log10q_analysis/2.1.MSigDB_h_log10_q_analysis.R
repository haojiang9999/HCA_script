#### 2.1.MSigDB_h_log10_q_analysis.R
COAD.tb.h.log10.q.all
library(reshape2)
library(ggplot2)
h.log10.q.all.m <- melt(COAD.tb.h.log10.q.all)
ggplot(h.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

rownames(COAD.tb.h.log10.q.all)











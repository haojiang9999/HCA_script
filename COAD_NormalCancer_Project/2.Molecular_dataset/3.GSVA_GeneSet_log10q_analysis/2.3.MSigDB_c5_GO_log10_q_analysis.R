#### 2.3.MSigDB_c5_GO_log10_q_analysis.R
c5.Names<- rownames(COAD.tb.c5.log10.q.all)
MSigDB.v7.0.sub.collection.Names <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/MSigDB/MSigDB.v7.0.sub.collection.Names.rds")
c5.names.sub <- MSigDB.v7.0.sub.collection.Names$c5.names
### BP: GO biological process
c5.log10.q.bp <- COAD.tb.c5.log10.q.all[c5.Names %in% c5.names.sub$c5.bp.names,]
dim(c5.log10.q.bp)
library(reshape2)
library(ggplot2)
c5.log10.q.bp.m <- melt(c5.log10.q.bp)
ggplot(c5.log10.q.bp.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "BP: GO biological process")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

#### 






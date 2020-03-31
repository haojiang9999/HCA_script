#### 2.2.MSigDB_c2_log10_q_analysis.R
c2.Names<- rownames(COAD.tb.c2.log10.q.all)
MSigDB.v7.0.sub.collection.Names <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/MSigDB/MSigDB.v7.0.sub.collection.Names.rds")
c2.names.sub <- MSigDB.v7.0.sub.collection.Names$c2.names
### CGP: chemical and genetic perturbations subset ###
c2.log10.q.cgp <- COAD.tb.c2.log10.q.all[c2.Names %in% c2.names.sub$c2.cgp.names,]
library(reshape2)
library(ggplot2)
c2.log10.q.cgp.m <- melt(c2.log10.q.cgp)
ggplot(c2.log10.q.cgp.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "CGP: chemical and genetic perturbations subset")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
rownames(c2.log10.q.cgp)[1:100]
#paste(rownames(c2.log10.q.cgp)[1:100],collapse = " AND ")
### KEGG subset ###
c2.log10.q.kegg <- COAD.tb.c2.log10.q.all[c2.Names %in% c2.names.sub$c2.cp.kegg.names,]
library(reshape2)
library(ggplot2)
c2.log10.q.kegg.m <- melt(c2.log10.q.kegg)
ggplot(c2.log10.q.kegg.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "KEGG pathway")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
rownames(c2.log10.q.kegg)[1:100]
## PATHWAY only
c2.log10.q.kegg.PATHWAY <- c2.log10.q.kegg[grep("PATHWAY",rownames(c2.log10.q.kegg)),]
c2.log10.q.kegg.PATHWAY.m <- melt(c2.log10.q.kegg.PATHWAY)
ggplot(c2.log10.q.kegg.PATHWAY.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "KEGG pathway")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
#### Metabolism ####  
c2.log10.q.kegg.METABOLISM <- c2.log10.q.kegg[grep("METABOLISM",rownames(c2.log10.q.kegg)),]
c2.log10.q.kegg.METABOLISM.m <- melt(c2.log10.q.kegg.METABOLISM)
ggplot(c2.log10.q.kegg.METABOLISM.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "KEGG METABOLISM")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()


### BIOCARTA subset ###
c2.log10.q.biocarta <- COAD.tb.c2.log10.q.all[c2.Names %in% c2.names.sub$c2.cp.biocarta.names,]
library(reshape2)
library(ggplot2)
c2.log10.q.biocarta.m <- melt(c2.log10.q.biocarta)
ggplot(c2.log10.q.biocarta.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + labs(title = "BIOCARTA pathwat")+
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
#### Activation ####
#c2.log10.q.biocarta.ACTIVATION <- c2.log10.q.biocarta[grep("ACTIVATION",rownames(c2.log10.q.biocarta)),]
#c2.log10.q.biocarta.ACTIVATION.m <- melt(c2.log10.q.biocarta.ACTIVATION)
#ggplot(c2.log10.q.biocarta.ACTIVATION.m, aes(Var2, Var1)) +
#  geom_tile(aes(fill = value),colour = "black") + labs(title = "biocarta.ACTIVATION")+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
#  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()



### REACTOME subset ###
c2.log10.q.reactome <- COAD.tb.c2.log10.q.all[c2.Names %in% c2.names.sub$c2.cp.reactome.names,]
dim(c2.log10.q.reactome)
library(reshape2)
library(ggplot2)
c2.log10.q.reactome.m <- melt(c2.log10.q.reactome)
ggplot(c2.log10.q.reactome.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + labs(title = "REACTOME pathwat")+
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

### PID subset ###
c2.log10.q.pid <- COAD.tb.c2.log10.q.all[c2.Names %in% c2.names.sub$c2.cp.pid.names,]
library(reshape2)
library(ggplot2)
c2.log10.q.pid.m <- melt(c2.log10.q.pid)
ggplot(c2.log10.q.pid.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + labs(title = "PID pathwat")+
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

#### Metabolism ####  
COAD.tb.c2.log10.q.all[grep("METABOLISM",c2.Names),]



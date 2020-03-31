#### 3.Gene_set_enrichment_analysis_Summary.R
### 1.Signatures
COAD.tb.h.log10.q.all <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/2.Molecular_dataset/2.Gene_set_analysis/COAD.tb.h.log10.q.all.rds")
## Names selection
h.names <- rownames(COAD.tb.h.log10.q.all)
h.selected <- c("HALLMARK_ALLOGRAFT_REJECTION",
                "HALLMARK_MYC_TARGETS_V2",
                "HALLMARK_KRAS_SIGNALING_UP",
                "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                "HALLMARK_IL2_STAT5_SIGNALING",
                "HALLMARK_COMPLEMENT",
                "EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_ANGIOGENESIS",
                "HALLMARK_APOPTOSIS",
                "HALLMARK_HYPOXIA",
                "GO_RESPONSE_TO_OXIDATIVE_STRESS",
                "SOUCEK_MYC_TARGETS")
library(reshape2)
library(ggplot2)
COAD.tb.h.log10.q.all.selected <- COAD.tb.h.log10.q.all[h.names %in% h.selected,1:3]
h.log10.q.all.m <- melt(COAD.tb.h.log10.q.all[h.names %in% h.selected,1:3])
ggplot(h.log10.q.all.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black",lwd = 1) + theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Hall markers seleted")+
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

write.csv(COAD.tb.h.log10.q.all.selected, file = "COAD.tb.h.log10.q.all.selected.csv")

### 2.Canonical KEGG pathway
## Loading
COAD.tb.c2.log10.q.all <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/2.Molecular_dataset/2.Gene_set_analysis/COAD.tb.c2.log10.q.all.rds")
###
c2.Names<- rownames(COAD.tb.c2.log10.q.all)
MSigDB.v7.0.sub.collection.Names <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/MSigDB/MSigDB.v7.0.sub.collection.Names.rds")
c2.names.sub <- MSigDB.v7.0.sub.collection.Names$c2.names
## KEGG pathway
c2.log10.q.kegg <- COAD.tb.c2.log10.q.all[c2.Names %in% c2.names.sub$c2.cp.kegg.names,]
library(reshape2)
library(ggplot2)

## PATHWAY only
c2.log10.q.kegg.PATHWAY <- c2.log10.q.kegg[grep("PATHWAY",rownames(c2.log10.q.kegg)),]
c2.log10.q.kegg.PATHWAY.m <- melt(c2.log10.q.kegg.PATHWAY)
ggplot(c2.log10.q.kegg.PATHWAY.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "KEGG pathway")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
rownames(c2.log10.q.kegg.PATHWAY)
## PATHWAY selection
KEGG.PATHWAY.selected <- c("KEGG_JAK_STAT_SIGNALING_PATHWAY",
                           "KEGG_MAPK_SIGNALING_PATHWAY",
                           "KEGG_CALCIUM_SIGNALING_PATHWAY",
                           "KEGG_VEGF_SIGNALING_PATHWAY",
                           "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                           "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
                           "KEGG_MTOR_SIGNALING_PATHWAY",
                           "KEGG_WNT_SIGNALING_PATHWAY",
                           "KEGG_NOTCH_SIGNALING_PATHWAY",
                           "KEGG_P53_SIGNALING_PATHWAY",
                           "PID_PI3KCI_PATHWAY",
                           "REACTOME_MTOR_SIGNALLING",
                           "GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
                           "BIOCARTA_AKT_PATHWAY",
                           "PID_HIF1_TFPATHWAY ")
c2.log10.q.kegg.PATHWAY.selected<- c2.log10.q.kegg.PATHWAY[rownames(c2.log10.q.kegg.PATHWAY) %in% KEGG.PATHWAY.selected,1:3]
c2.log10.q.kegg.PATHWAY.m <- melt(c2.log10.q.kegg.PATHWAY[rownames(c2.log10.q.kegg.PATHWAY) %in% KEGG.PATHWAY.selected,1:3])
ggplot(c2.log10.q.kegg.PATHWAY.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(title = "KEGG pathway seleted")+
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
write.csv(c2.log10.q.kegg.PATHWAY.selected,file = "c2.log10.q.kegg.PATHWAY.selected.csv")


#### Metabolism ####  
c2.log10.q.kegg.METABOLISM <- c2.log10.q.kegg[grep("METABOLISM",rownames(c2.log10.q.kegg)),]
c2.log10.q.kegg.METABOLISM.m <- melt(c2.log10.q.kegg.METABOLISM)
ggplot(c2.log10.q.kegg.METABOLISM.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + labs(title = "KEGG METABOLISM")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()
rownames(c2.log10.q.kegg.METABOLISM)
## Metabolism selection
KEGG.METABOLISM.selected <- c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                              "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM",
                              "KEGG_STARCH_AND_SUCROSE_METABOLISM",
                              "KEGG_GALACTOSE_METABOLISM",
                              "KEGG_GLUTATHIONE_METABOLISM",
                              "KEGG_NITROGEN_METABOLISM",
                              "KEGG_TYROSINE_METABOLISM",
                              "KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
                              "KEGG_FATTY_ACID_METABOLISM",
                              "KEGG_LINOLEIC_ACID_METABOLISM",
                              "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
                              "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
                              "KEGG_OXIDATIVE_PHOSPHORYLATION",
                              "REACTOME_BETA_OXIDATION_OF_VERY_LONG_CHAIN_FATTY_ACIDS",
                              "REACTOME_FATTY_ACIDS",
                              "GO_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_OXIDASE",
                              "KEGG_PYRUVATE_METABOLISM"
                              )

c2.log10.q.kegg.METABOLISM.selected <- c2.log10.q.kegg[rownames(c2.log10.q.kegg) %in% KEGG.METABOLISM.selected,1:3]
c2.log10.q.kegg.METABOLISM.selected.m <- melt(c2.log10.q.kegg.METABOLISM.selected)
ggplot(c2.log10.q.kegg.METABOLISM.selected.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "black") + theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "KEGG metabolism seleted")+
  scale_fill_gradient2(low = "mediumblue", high = "red2", mid = "white",midpoint = 0) + coord_fixed()

write.csv(c2.log10.q.kegg.METABOLISM.selected, file = "c2.log10.q.kegg.METABOLISM.selected.csv")

### Immune selection




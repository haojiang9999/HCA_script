#### 3.1.Pathway_selection.R
### 1.Loading
COAD.tb.h.log10.q.all <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/2.Molecular_dataset/2.Gene_set_analysis/COAD.tb.h.log10.q.all.rds")
COAD.tb.c2.log10.q.all <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/2.Molecular_dataset/2.Gene_set_analysis/COAD.tb.c2.log10.q.all.rds")
COAD.tb.c5.log10.q.all <- readRDS("/data8t_4/JH/MyJobs/COAD_NormalCancer_Project/2.Molecular_dataset/2.Gene_set_analysis/COAD.tb.c5.log10.q.all.rds")
COAD.tb.h.c2.c5.log10.q <- rbind(COAD.tb.h.log10.q.all, COAD.tb.c2.log10.q.all, COAD.tb.c5.log10.q.all)

### Metabolism ###
anaplerotic ANAPLER # none
lactate LACTA
glutaminolysis GLUTAMINO
glycolysis GLYCOLYSIS
oxidative OXIDATIVE
glutaminolysis GLUTAMINO # none
fatty acid synthesis FATTY
pyruvate PYRUVATE
redox REDOX # none in c2 less in c5
lipid synthesis LIPID_SYN # none
synthesis SYNTHESIS
## h c2 c5
COAD.tb.h.c2.c5.log10.q[grep("OXIDATIVE",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("FATTY",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("PYRUVATE",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("PROTEIN",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("GLUTAMINE",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]

head(COAD.tb.h.c2.c5.log10.q[grep("PROTEIN",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F])
### PATHWATS ###
PI3K
COAD.tb.h.c2.c5.log10.q[grep("PI3K",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
MTOR
COAD.tb.h.c2.c5.log10.q[grep("MTOR",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
AKT
COAD.tb.h.c2.c5.log10.q[grep("AKT",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("KRAS",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("NF_KB",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("TGF_BETA",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
### signatures
MYC
COAD.tb.h.c2.c5.log10.q[grep("MYC",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
RAS
COAD.tb.h.c2.c5.log10.q[grep("RAS",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
HIF
COAD.tb.h.c2.c5.log10.q[grep("HIF",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("EMT",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("MATRIX",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("AUTO",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
### Immune states
## Inflammation
COAD.tb.h.c2.c5.log10.q[grep("IFN",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("MYD88",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("IL6",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("IL17",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]
COAD.tb.h.c2.c5.log10.q[grep("BETA",rownames(COAD.tb.h.c2.c5.log10.q)),, drop = F]

#### Pipline for Organoid analysis
# Orgnoid_expression.uniGene.TPM
### Step1.expression data transformation
log2.Organoid_expression <- log2(Orgnoid_expression.uniGene.TPM +1)
### Step2.Distance calculation
##### Distance calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
log2.Organoid_expression.list <- list(log2.Organoid_expression = log2.Organoid_expression)

Cor.Res.CV8000.Org <- refCorMerge(log2.Organoid_expression.list, log2.scReference.list.CV.8000)
Cor.Res.CV4000.Org <- refCorMerge(log2.Organoid_expression.list, log2.scReference.list.CV.4000)
Cor.Res.CV2000.Org <- refCorMerge(log2.Organoid_expression.list, log2.scReference.list.CV.2000)
Cor.Res.CV1500.Org <- refCorMerge(log2.Organoid_expression.list, log2.scReference.list.CV.1500)
Cor.Res.CV1000.Org <- refCorMerge(log2.Organoid_expression.list, log2.scReference.list.CV.1000)
### Pheno info merge
Pheno.merged.Organoid <- Cor.Res.CV1000.Org$Pheno.merged
drug.res <- Orgnoid.colon.AUC.IC50$drug.res
colnames(drug.res)
rownames(Pheno.merged.Organoid)
Pheno.merged.Organoid <- cbind(Pheno.merged.Organoid, t(drug.res))
###
### Step3.Dimension reduction analysis
source("/data8t_4/JH/MyJobs/1_R_script/R_Plot/DimeReduPlot.R")
DimeReduPlot(mx = scale(Cor.Res.CV8000.Org$Cor.merged), color = Pheno.merged.Organoid$X5_FU_AUC, 
             tiltle = "Cor.Res.CV8000.Org_X5_FU_AUC", print = T, perplexity = 10)

DimeReduPlot(mx = scale(Cor.Res.CV8000.Org$Cor.merged), color = Pheno.merged.Organoid$Nutlin3a_AUC, 
             tiltle = "Cor.Res.CV8000.Org_X5_FU_AUC", print = T, perplexity = 10)

DimeReduPlot(mx = scale(Cor.Res.CV2000.Org$Cor.merged), color = Pheno.merged.Organoid$Afatinib_IC50, 
             tiltle = "Cor.Res.CV2000.Org_X5_FU_AUC", print = T, perplexity = 10)

DimeReduPlot(mx = scale(Cor.Res.CV1500.Org$Cor.merged), color = Pheno.merged.Organoid$Doxorubicin_IC50, 
             tiltle = "Cor.Res.CV1500.Org_X5_FU_AUC", print = T, perplexity = 10)
## Heatmap

ClustHeatmap( Cor.Res.CV1500.Org$Cor.merged, Pheno.merged= Pheno.merged.Organoid, title = "Cor.Res.CV1500.Org",scale = c("column"))
ClustHeatmap( Cor.Res.CV2000.Org$Cor.merged, Pheno.merged= Pheno.merged.Organoid, title = "Cor.Res.CV2000.Org",scale = c("column"))
ClustHeatmap( Cor.Res.CV4000.Org$Cor.merged, Pheno.merged= Pheno.merged.Organoid, title = "Cor.Res.CV4000.Org",scale = c("column"))


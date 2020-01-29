##### Distance calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
Cor.Res.CV8000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.8000)
Cor.Res.CV4000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.4000)
Cor.Res.CV2000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.2000)
Cor.Res.CV1500 <- refCorMerge(Log2.expList, log2.scReference.list.CV.1500)
Cor.Res.CV1000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.1000)


### Need 3.1 PhenoType information

Log2.Cor.Res.list <- list(Cor.Res.CV8000 = Cor.Res.CV8000,
                          Cor.Res.CV4000 = Cor.Res.CV4000,
                          Cor.Res.CV2000 = Cor.Res.CV2000,
                          Cor.Res.CV1500 = Cor.Res.CV1500,
                          Cor.Res.CV1000 = Cor.Res.CV1000,
                          Pheno.merged = Pheno.merged)
saveRDS(Log2.Cor.Res.list, file = "Log2.Cor.Res.list.dataset.rds")

### Read saved data
Log2.Cor.Res.list <- readRDS("Log2.Cor.Res.list.dataset.rds")
Cor.Res.CV2000 <- Log2.Cor.Res.list$Cor.Res.CV2000
Cor.Res.CV1500 <- Log2.Cor.Res.list$Cor.Res.CV1500
Pheno.merged <- Log2.Cor.Res.list$Pheno.merged

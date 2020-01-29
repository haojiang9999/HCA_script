### Distance calculation
##### Distance calculation 
source("/data8t_4/JH/MyJobs/1_R_script/NormalCancer/refCorMerge.R")
### 
##### Transform COAD data
summary(colSums(COAD_tpm_symbol)) ### So its not log2 transformed
# Expression data transformation Log(2+0.001)
Log2.expList <- list(COAD_tpm_symbol = log2(COAD_tpm_symbol+1))
Cor.Res.CV8000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.8000)
Cor.Res.CV4000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.4000)
Cor.Res.CV2000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.2000)
Cor.Res.CV1500 <- refCorMerge(Log2.expList, log2.scReference.list.CV.1500)
Cor.Res.CV1000 <- refCorMerge(Log2.expList, log2.scReference.list.CV.1000)

### Check out
summary(Cor.Res.CV8000$Cor.merged)
summary(Cor.Res.CV1500$Cor.merged)

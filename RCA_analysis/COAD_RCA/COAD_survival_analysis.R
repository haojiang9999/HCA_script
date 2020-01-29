#### COAD survival analysis ####
## Step1:build data frame
# cluster results
COAD.cor.cluster <- clusterResault[[2]]
## Step2: read pateint survival data
COAD_UCSC_Toil_tpm_dataset <- readRDS("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/COAD_UCSC_Toil_tpm_dataset.rds")
COAD.pheno.exp <- COAD_UCSC_Toil_tpm_dataset$COAD.pheno.exp
## sample type ###
sampleName <- rownames(COAD.pheno.exp)
samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
table(samplType)
tumorIndex <- samplType == "01"

## Step3: build data frame
sampleID <- rownames(COAD.cor.cluster)
surTime <- COAD.pheno.exp[sampleID,c("OS","OS.time")]
COAD.sur.df <- cbind(COAD.cor.cluster,surTime)

### survival analysis
library(survminer)
require("survival")
fit <- survfit(Surv(OS.time, OS) ~ dynamicColors, data = COAD.sur.df)
# Drawing curves
ggsurvplot(fit)

#### remove slusters have little samples
biggroup <- table(COAD.sur.df$dynamicColors) > 20
colornames <- names(biggroup[biggroup == T])
COAD.sur.df.test <- COAD.sur.df[COAD.sur.df$dynamicColors %in% colornames,]

fit <- survfit(Surv(OS.time, OS) ~ dynamicColors, data = COAD.sur.df.test)
# Drawing curves
ggsurvplot(fit)















#### 3.Read_COADREAD_Nature_2012.maf
## Download : https://www.synapse.org/#!Synapse:syn2812925
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")
library(maftools)
COAD.maf <- "COADREAD_Nature_2012.maf"
COAD <- read.maf(maf = COAD.maf)
COAD
#Shows gene summary.
getGeneSummary(COAD)
#Shows sample summry.
getSampleSummary(COAD)
#Shows all fields in MAF
getFields(COAD)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = COAD, basename = 'COAD')9

#oncoplot for top ten mutated genes.
oncoplot(maf = COAD, top = 10)

oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))

# read gene annotation file
gtfPath <- "/stor/jianghao/GEO/Normal_tissues/GSE103239_Tang_digestive_tract/GSE103154_adult_tissues/UCSC_hg19_gene_annotation.gtf.gz"
UCSC.hg19.gtf <- read.table(gtfPath)
head(UCSC.hg19.gtf)
table(UCSC.hg19.gtf$V2)
UCSC.hg19.gtf$V10
#
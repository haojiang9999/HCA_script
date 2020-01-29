#### Generate Xena TCGA methylation dataset
# https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#### 1.read RNA-seq expression data sets ####
filePath <- "/stor/jianghao/Xena/TCGA_Pan_Cancer/DNA_methylation/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz"
#### Data (file names: *.rsem_genes.results) are downloaded, 
#    tpm values are extracted, log2(x+0.001) transformed
PANCAN_HumanMethylation450.betaValue <- read.delim(filePath, header = T)







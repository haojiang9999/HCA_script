#### 2.Geneset_list_construction.R
### Read data from https://www.synapse.org/#!Synapse:syn2623706/wiki/67246
# The consensus molecular subtypes of colorectal cancer
# Consensus Molecular Subtypes(CMS)
genesets.cms <- read.delim("/stor/jianghao/Synapse/geneset/genesets.gmt",header = F)
genesets.cms.names <- as.character(genesets.cms$V1)
genesets.cms.list <- split(genesets.cms[,-1], seq(nrow(genesets.cms[,-1])))
### Remove "" empty elements in the lists
genesets.cms.list<- lapply(genesets.cms.list, function(x){
  x[] <- lapply(x, as.character)
  geneNames <- as.character(x)
  geneNames <- geneNames[geneNames != ""]
  return(geneNames)
})
#### Compare gmt genes in the TCGA expression table
allgenes <- unlist(genesets.cms.list, use.names=FALSE) 
allgenes <- unique(allgenes)
TCGA.genes <- rownames(COAD.RSEM.gene.norm.count.Log2)
TCGA.genes2 <- rownames(COAD.RSEM.gene.tpm.symbol)
table(allgenes %in% TCGA.genes)
table(allgenes %in% TCGA.genes2)
allgenes[!allgenes %in% TCGA.genes]
grep("TXNDC7",TCGA.genes )
#### Shouled I further to find the unmaped 420 allgenes


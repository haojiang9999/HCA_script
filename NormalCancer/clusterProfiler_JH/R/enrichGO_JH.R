##' Modified GO Enrichment Analysis from clusterprofiler enrichGO.
##' 
##' 
##' @description 
##' This modified fucntion through pre-load GO_DATA("CC","MF","BP) 
##' into the enviroment to save time.
##' Given a vector of genes, this function will return the enrichment GO
##' categories after FDR control.
##' @param gene a vector of entrez gene id.
##' @param OrgDb OrgDb
##' @param keyType keytype of input gene
##' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param readable whether mapping gene ID to gene Name
##' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
##' @return An \code{enrichResult} instance.
##' @importClassesFrom DOSE enrichResult
##' @importFrom DOSE setReadable
##' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @examples
##' \dontrun{
##'   data(geneList, package = "DOSE")
##' 	de <- names(geneList)[1:100]
##' 	yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
##' 	head(yy)
##' }
##' # have to keep rds in some local folder
##' # "/stor/jianghao/Database/GO_DATA"
##' GO_DATA_BP <- readRDS("./GO_DATA/GO_DATA_BP.rds")
##' GO_DATA_CC <- readRDS("./GO_DATA/GO_DATA_CC.rds")
##' GO_DATA_MF <- readRDS("./GO_DATA/GO_DATA_MF.rds")
##'
##' ego_CC <- enrichGO_JH(gene          = gene,
##'                      universe      = names(geneList),
##'                      OrgDb         = org.Hs.eg.db,
##'                      ont           = "CC", # chose ontology terms
##'                      pAdjustMethod = "BH", 
##'                      pvalueCutoff  = 0.01,
##'                      qvalueCutoff  = 0.05,
##'                      readable      = TRUE)
##' head(ego_BP)
##'                      
##' ego_CC <- enrichGO_JH(gene          = gene,
##' universe      = names(geneList),
##' OrgDb         = org.Hs.eg.db,
##' ont           = "BP", # chose ontology terms
##' pAdjustMethod = "BH", 
##' pvalueCutoff  = 0.01,
##' qvalueCutoff  = 0.05,
##' readable      = TRUE)                     
##'                      
##'                   
enrichGO_JH <- function (gene, OrgDb, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, 
          pAdjustMethod = "BH", universe, qvalueCutoff = 0.2, minGSSize = 10, 
          maxGSSize = 500, readable = FALSE, pool = FALSE) 
{
  #ont %<>% toupper
  ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
  if (ont == "CC"){
    GO_DATA <- GO_DATA_CC
  }else if (ont == "BP"){
    GO_DATA <- GO_DATA_BP
  }else if (ont == "MF"){
    GO_DATA <- GO_DATA_MF
  }
  
  if (missing(universe)) 
    universe <- NULL
  if (ont == "ALL" && !pool) {
    lres <- lapply(c("BP", "CC", "MF"), function(ont) suppressMessages(enrichGO(gene, 
                                                                                OrgDb, keyType, ont, pvalueCutoff, pAdjustMethod, 
                                                                                universe, qvalueCutoff, minGSSize, maxGSSize)))
    lres <- lres[!sapply(lres, is.null)]
    if (length(lres) == 0) 
      return(NULL)
    df <- do.call("rbind", lapply(lres, as.data.frame))
    geneSets <- lres[[1]]@geneSets
    if (length(lres) > 1) {
      for (i in 2:length(lres)) {
        geneSets <- append(geneSets, lres[[i]]@geneSets)
      }
    }
    res <- lres[[1]]
    res@result <- df
    res@geneSets <- geneSets
  }
  else {
    res <- enricher_internal(gene, pvalueCutoff = pvalueCutoff, 
                             pAdjustMethod = pAdjustMethod, universe = universe, 
                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, 
                             maxGSSize = maxGSSize, USER_DATA = GO_DATA)
    if (is.null(res)) 
      return(res)
  }
  res@keytype <- keyType
  res@organism <- get_organism(OrgDb)
  if (readable) {
    res <- setReadable(res, OrgDb)
  }
  res@ontology <- ont
  if (ont == "ALL") {
    res <- add_GO_Ontology(res, GO_DATA)
  }
  return(res)
}
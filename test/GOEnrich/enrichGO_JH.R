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
<bytecode: 0x62ed10a8>
  <environment: namespace:clusterProfiler>
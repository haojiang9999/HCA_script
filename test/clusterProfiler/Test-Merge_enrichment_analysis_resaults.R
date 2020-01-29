### test EnAnno
test7 <- EnAnno(inputID)
EnAnno.time1 <- system.time(test7 <- EnAnno(inputID))
EnAnno.time2 <- system.time(test7 <- EnAnno(inputID))

library(clusterProfiler)
test8.time <- system.time(GOEnrich(inputID))
enrichGO
OrgDb <- "org.Hs.eg.db"
OrgDb<- GOSemSim::load_OrgDb(OrgDb)
kt <- AnnotationDbi::keytypes(OrgDb)
if (! keytype %in% kt) {
  stop("keytype is not supported...")
}

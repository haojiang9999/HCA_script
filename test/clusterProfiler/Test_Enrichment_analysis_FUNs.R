
# test WikiEnrich 
inputID <- geneSymbol2GeneID(input)
inputID <- inputID$ENTREZID
wpgmtfile <- "/stor/jianghao/Database/WikiPathways/wikipathways-20190810-gmt-Homo_sapiens.gmt"

test6 <- WikiEnrich(inputID, wpgmtfile)### you background genes

head(test6@result)
DT::datatable(as.data.frame(test6))

#test CellMarkerEnrich
cellMarkerFile <- "/stor/jianghao/Database/Cell_Marker/Human_cell_markers.txt"

test5 <- CellMarkerEnrich( inputID, cellMarkerFile)
DT::datatable(as.data.frame(test5))  
test2@result


#GOEnrich
test<- GOEnrich(inputID)
head(test$egoCC)
head(test$egoMF)
head(test$egoBP)
head(test)

## KEGGEnrich

test2 <- KEGGEnrich(inputID)
head(test2$kk)
head(test2$mkk)


## ReactomeEnrich
test3<- ReactomeEnrich(inputID)
head(test3)

##MSigDBEnrich
test4<- MSigDBEnrich(inputID)
head(test4$C3_motif_gene_sets)
head(test4$C6_oncogenic_signaturess)
head(test4$C7_immunologic_signatures)
head(test4$H_hallmark_gene_sets)

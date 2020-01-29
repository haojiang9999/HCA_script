#### Step1.Generate enrichment resaults ####
# 1.WikiEnrich
wpgmtfile <- "/stor/jianghao/Database/WikiPathways/wikipathways-20190810-gmt-Homo_sapiens.gmt"
wikiEn <- WikiEnrich(inputID, wpgmtfile)
# 2.CellMarkerEnrich
cellMarkerFile <- "/stor/jianghao/Database/Cell_Marker/Human_cell_markers.txt"
cellEn <- CellMarkerEnrich( inputID, cellMarkerFile)
# 3.MSigDBEnrich
msigdbFilePath <- "/stor/jianghao/Database/MSigDB/"
MSEn <- MSigDBEnrich(inputID,msigdbFilePath = msigdbFilePath)
MSigDBEnrich.time <- system.time(MSEn <- MSigDBEnrich(inputID,msigdbFilePath = msigdbFilePath))
# 4.GOEnrich
GOEn <- GOEnrich(inputID)
# 5.KEGGEnrich
KEGGEn <- KEGGEnrich(inputID)
# 6.ReactomeEnrich
ReEn <- ReactomeEnrich(inputID)
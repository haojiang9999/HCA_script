#### merge enrichment analysis resaults ####
EnAnno <- function(inputID, background, TopN = 10){
  if (!require("tidyr")) 
    BiocManager::install("tidyr")
  #### Step1.Generate enrichment resaults ####
  # 1.WikiEnrich
  #wpgmtfile <- "/stor/jianghao/Database/WikiPathways/wikipathways-20190810-gmt-Homo_sapiens.gmt"
  wikiEn <- WikiEnrich(inputID, wpgmtfile, background = background)
  # 2.CellMarkerEnrich
  #cellMarkerFile <- "/stor/jianghao/Database/Cell_Marker/Human_cell_markers.txt"
  cellEn <- CellMarkerEnrich( inputID, cellMarkerFile, background = background)
  # 3.MSigDBEnrich
  msigdbFilePath <- "/stor/jianghao/Database/MSigDB/"
  MSEn <- MSigDBEnrich(inputID,msigdbFilePath = msigdbFilePath, background = background)
  # 4.GOEnrich
  GOEn <- GOEnrich(inputID, background = background)
  # 5.KEGGEnrich
  KEGGEn <- KEGGEnrich(inputID, background = background)
  # 6.ReactomeEnrich
  ReEn <- ReactomeEnrich(inputID, background = background)
  
  #### Step2.select top enrich resaults ####
  ## ordered by p-value
  EnrichRes <- list()
  #TopN <- 10
  ### Enrichment top selection
  # 1.wiki pathway
  EnrichRes$wikiEnTop <- head(wikiEn, n = TopN)
  # 2.cellType enrich
  EnrichRes$cellEnTop <- head(cellEn, n = TopN)
  # 3.msigdb enrich
  #looks no meaning
  #EnrichRes$MSEn_C1_positional_gene_setsTop <- head(MSEn$C1_positional_gene_sets, n = TopN)
  EnrichRes$MSEn_C2_curated_gene_setsTop <- head(MSEn$C2_curated_gene_sets, n = TopN)
  EnrichRes$MSEn_C3_motif_gene_setsTop <- head(MSEn$C3_motif_gene_sets, n = TopN)
  EnrichRes$MSEn_C4_computational_gene_setsTop <- head(MSEn$C4_computational_gene_sets, n = TopN)
  # have another GO analysis
  #EnrichRes$MSEn_C5_GO_gene_setsTop <- head(MSEn$C5_GO_gene_sets, n = TopN)
  EnrichRes$MSEn_C6_oncogenic_signaturesTop <- head(MSEn$C6_oncogenic_signatures, n = TopN)
  EnrichRes$MSEn_C7_immunologic_signaturesTop <- head(MSEn$C7_immunologic_signatures, n = TopN)
  EnrichRes$MSEn_H_hallmark_gene_setsTop <- head(MSEn$H_hallmark_gene_sets, n = TopN)
  # 4.GOEnrich
  EnrichRes$GOEnTop_CC <- head(GOEn$egoCC, n = TopN)
  EnrichRes$GOEnTop_MF <- head(GOEn$egoMF, n = TopN)
  EnrichRes$GOEnTop_BP <- head(GOEn$egoBP, n = TopN)
  # 5.KEGGEn
  #EnrichRes$KEGGEnTop_kk <- head(KEGGEn$kk, n = TopN)
  #EnrichRes$KEGGEnTop_mkk <- head(KEGGEn$mkk, n = TopN)
  # 6.ReactomeEnrich
  EnrichRes$ReEnTop <- head(ReEn, n = TopN)
  return(EnrichRes)
  
} 

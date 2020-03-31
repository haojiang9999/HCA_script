#### 3.Sub_collections_annatation.R
library(GSEABase)
#### c2 sub collection names #####
c2.cgp.v7.0.symbols <- getGmt("c2_sub_collection/c2.cgp.v7.0.symbols.gmt")
c2.cp.biocarta.v7.0.symbols <- getGmt("c2_sub_collection/c2.cp.biocarta.v7.0.symbols.gmt")
c2.cp.kegg.v7.0.symbols <- getGmt("c2_sub_collection/c2.cp.kegg.v7.0.symbols.gmt")
c2.cp.pid.v7.0.symbols <- getGmt("c2_sub_collection/c2.cp.pid.v7.0.symbols.gmt")
c2.cp.reactome.v7.0.symbols <- getGmt("c2_sub_collection/c2.cp.reactome.v7.0.symbols.gmt")
c2.cp.v7.0.symbols <- getGmt("c2_sub_collection/c2.cp.v7.0.symbols.gmt")
### Extract names
c2.cgp.names<- names(c2.cgp.v7.0.symbols)
c2.cp.biocarta.names<- names(c2.cp.biocarta.v7.0.symbols)
c2.cp.kegg.names<- names(c2.cp.kegg.v7.0.symbols)
c2.cp.pid.names<- names(c2.cp.pid.v7.0.symbols)
c2.cp.reactome.names<- names(c2.cp.reactome.v7.0.symbols)
c2.cp.names<- names(c2.cp.v7.0.symbols)

c2.names <- list(c2.cgp.names = c2.cgp.names, c2.cp.biocarta.names = c2.cp.biocarta.names,
     c2.cp.kegg.names = c2.cp.kegg.names,c2.cp.pid.names = c2.cp.pid.names,
     c2.cp.reactome.names = c2.cp.reactome.names, c2.cp.names = c2.cp.names,
     CGP = "CGP: chemical and genetic perturbations",
     CP = "CP: Canonical pathways")

#### c3 sub collection names #####
c3.mir.v7.0.symbols <- getGmt("c3_sub_collection/c3.mir.v7.0.symbols.gmt")
c3.tft.v7.0.symbols <- getGmt("c3_sub_collection/c3.tft.v7.0.symbols.gmt")
### Extract names
c3.mir.names <- names(c3.mir.v7.0.symbols)
c3.tft.names <- names(c3.tft.v7.0.symbols)

c3.names <- list(c3.mir.names = c3.mir.names, c3.tft.names = c3.tft.names, 
                 MIR = "MIR: microRNA targets",
                 TFT = "TFT: transcription factor targets")

#### c4 sub collection names #####
c4.cgn.v7.0.symbols <- getGmt("c4_sub_collection/c4.cgn.v7.0.symbols.gmt")
c4.cm.v7.0.symbols <- getGmt("c4_sub_collection/c4.cm.v7.0.symbols.gmt")

c4.cgn.names<- names(c4.cgn.v7.0.symbols)
c4.cm.names<- names(c4.cm.v7.0.symbols)

c4.names <- list(c4.cgn.names = c4.cgn.names, c4.cm.names = c4.cm.names, 
                 CGN = "CGN: cancer gene neighborhoods",
                 CM = "CM: cancer modules")

#### c5 sub collection names #####
c5.bp.v7.0.symbols <- getGmt("c5_sub_collection/c5.bp.v7.0.symbols.gmt")
c5.cc.v7.0.symbols <- getGmt("c5_sub_collection/c5.cc.v7.0.symbols.gmt")
c5.mf.v7.0.symbols <- getGmt("c5_sub_collection/c5.mf.v7.0.symbols.gmt")
#\
c5.bp.names <- names(c5.bp.v7.0.symbols)
c5.cc.names <- names(c5.cc.v7.0.symbols)
c5.mf.names <- names(c5.mf.v7.0.symbols)

c5.names <- list(c5.bp.names = c5.bp.names, c5.cc.names = c5.cc.names,
                 c5.mf.names = c5.mf.names, BP= "BP: GO biological process",
                 CC = "CC: GO cellular component", 
                 MF = "MF: GO molecular function")

MSigDB.v7.0.sub.collection.Names <- list(c2.names = c2.names, c3.names = c3.names, c4.names = c4.names,
                                         c5.names = c5.names)
saveRDS(MSigDB.v7.0.sub.collection.Names, file = "MSigDB.v7.0.sub.collection.Names.rds")

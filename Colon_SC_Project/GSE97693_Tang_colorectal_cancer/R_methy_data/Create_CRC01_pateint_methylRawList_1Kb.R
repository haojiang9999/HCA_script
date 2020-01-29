#CRC01 pateint
#Genomes were divided into 1Kb tiles and genome methylation level 
#were also measured by 1Kb
#creat first 10 samples in methylRawList
library(methylKit)
filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/GSE97693_Met"
file.CpG <- paste0(filePath,'/',sample.CRC01[1:2],".CpG.txt.gz")
objDB.CRC01 <- methRead(as.list(file.CpG),sample.id=as.list(CRC01.id[1:2]),
                        treatment = c(1,1),
                        mincov=1, # this was important set to 1 in single cell data 
                        assembly="hg9",header=F, context="CpG", resolution="base",
                        pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                      coverage.col=5,strand.col=4,freqC.col=8 ))



#read files as for loop frome Num 2 sample
for (i in 1:length(sample.CRC01)) {
  start_time <- Sys.time()
  filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/GSE97693_Met"
  file.CpG <- paste0(filePath,'/',sample.CRC01[i],".CpG.txt.gz")
  #  print(file.CpG)
  sampleId <- as.character(sample.CRC01[i])
  print(sampleId)
  # 1.create methylKit object
  obj <- methRead(file.CpG,sample.id=CRC01.id[i],
                  mincov=1, # this was important set to 1 in single cell data 
                  assembly="hg9",header=F, context="CpG", resolution="base",
                  pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                coverage.col=5,strand.col=4,freqC.col=8 ))
  # 2.Tiling windows analysis
  tile <- tileMethylCounts(obj,win.size=1000,step.size=1000)
  #print(tile)
  # 3.add a new methylRaw object to methylRawList
  objDB.CRC01[[i]] <- tile
  # 4.update the treatment vector in the methylRawList object
  #    the treatment vector designates the sample groups
  objDB.CRC01@treatment=c(objDB.CRC01@treatment,1)
  end_time <- Sys.time()
  print(end_time - start_time)
}
length(objDB.CRC01)
format(object.size(objDB.CRC01[[1]]), units = "MB")
objDB.CRC01
######### continue where left ############
for (i in 255:length(sample.CRC01)) {
  start_time <- Sys.time()
  filePath <- "/data8t_4/JH/scRNA_seq/GEO/Cancer/GSE97693_Tang_colorectal_cancer/GSE97693_Met"
  file.CpG <- paste0(filePath,'/',sample.CRC01[i],".CpG.txt.gz")
  #  print(file.CpG)
  sampleId <- as.character(sample.CRC01[i])
  print(sampleId)
  # 1.create methylKit object
  obj <- methRead(file.CpG,sample.id=CRC01.id[i],
                  mincov=1, # this was important set to 1 in single cell data 
                  assembly="hg9",header=F, context="CpG", resolution="base",
                  pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                coverage.col=5,strand.col=4,freqC.col=8 ))
  # 2.Tiling windows analysis
  tile <- tileMethylCounts(obj,win.size=1000,step.size=1000)
  #print(tile)
  # 3.add a new methylRaw object to methylRawList
  objDB.CRC01[[i]] <- tile
  # 4.update the treatment vector in the methylRawList object
  #    the treatment vector designates the sample groups
  objDB.CRC01@treatment=c(objDB.CRC01@treatment,1)
  end_time <- Sys.time()
  print(end_time - start_time)
}
objDB.CRC01[[408]]
saveRDS(objDB.CRC01, file = "MethyKit.objDB.CRC01.409cells.1Kb.rds")

rm(objDB.CRC01)


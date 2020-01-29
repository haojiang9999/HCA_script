#### test for TCGA.surv ####
TCGA.surv(cancerType = "COAD",geneName = DE_lncRNAs[4],
          print = T)
TCGA.surv()
system.time(TCGA.surv(cancerType = "COAD",geneName = "ENSG00000172965",
                      print = T))
#### TCGA cancer types ##### 
TCGA.files <- list.files("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/")
TCGA.files <- TCGA.files[grep("dataset.rds", TCGA.files)]
cancerTypes.short <- unlist(lapply(strsplit(TCGA.files,"_"), '[[', 1))
sapply(Gene_ID[1:3],FUN = function(x){
  print(x)
  #TCGA.surv(cancerType = "COAD",geneName = x,print = T)
  TCGA.expCompare(cancerType = "COAD",geneName = "ENSG00000172965",print = T)
})

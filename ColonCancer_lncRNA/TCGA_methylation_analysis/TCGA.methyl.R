#### Function TCGA.methyl.R
TCGA.methyl <- function(cancerType = "COAD",geneName = "C7orf40",print = TRUE){
  #### 1.Reaad TCGA dataset
  filePath <- paste0("/stor/jianghao/Myjobs/R/Read_data_set/TCGA_DNA_mythelation/",cancerType,"_PANCAN_Methylation450_dataset.rds")
  TCGA.Methylation450.dataset <- readRDS(filePath)
  #TCGA.Methylation450.dataset$COADMethylation450
  ## Extract Methylation data
  TCGA.Methylation450 <- TCGA.Methylation450.dataset[[paste0(cancerType,"Methylation450")]]
  ## Gene and cg#### number annotation data
  Methyl450_hg19_GPL16304 <- TCGA.Methylation450.dataset$Methyl450_hg19_GPL16304
  ## Extract gene correlated cg####s
  geneCGIndex <- c(grep(geneName, Methyl450_hg19_GPL16304$gene),
                   grep(geneName, Methyl450_hg19_GPL16304$X.id))
  cgNames <- as.character(Methyl450_hg19_GPL16304[geneCGIndex,]$X.id)
  if(identical(cgNames, character(0))){
    print(paste0(geneName," not exits in TCGA methylation dataset."))
  }else{
    ## Gene cg##### expression matrix
    TCGA.Methylation450.sub.gene <- TCGA.Methylation450[cgNames,]
    #### 2.Separete TCGA data into normal and cancer types
    sampleName <- colnames(TCGA.Methylation450.sub.gene)
    samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
    #table(samplType)
    ## "01" was primary tumor, "11' was normal tissue
    tumorIndex <- samplType == "01"
    normalIndex <- samplType == "11"
    ## Expression data for normal and cancer
    tumor.meth <- TCGA.Methylation450.sub.gene[, tumorIndex]
    normal.meth<- TCGA.Methylation450.sub.gene[, normalIndex]
    # Number of tumor and normal samples
    NumSam.tumor <- length(colnames(tumor.meth))
    NumSam.normal <- length(colnames(normal.meth))
    #### 3.Build ggplot dataframe list
    df <- rbind(t(tumor.meth),t(normal.meth))
    df <- data.frame(df, TissueTypes = c(rep("Tumor",NumSam.tumor),rep("Normal",NumSam.normal)))
    # check
    #df[df$TissueTypes == "Tumor",]
    #head(tumor.meth)
    #### 4.Plot
    #i=7
    res <- list()
    for(i in 1:length(cgNames)){
      # t-test
      p_value <- t.test(df[,cgNames[i]] ~ df$TissueTypes )
      # figure plot
      p <- ggplot(df, aes_string("TissueTypes",cgNames[i])) # In here using aes_string!!!
      res[[i]] <- p + geom_boxplot(fill = c("dodgerblue3","firebrick3")) +  
        labs(x= paste0("Normal.N = ",NumSam.normal," Tumor.N = ",NumSam.tumor),
             y = "DNA methylation Beta value",
             title =paste0(geneName," in ", cancerType," at ",cgNames[i]),
             subtitle = paste0("p-value = ", p_value$p.value)
        )
      
    }
    if(print){
      library(ggpubr)
      print(ggarrange(plotlist = res))
    }else{
      res.arr <- ggarrange(plotlist = res, ncol = 2)
      ggexport(res.arr, filename = paste0(geneName," in ",cancerType," methylation analysis.pdf"))
    } 
  }
  
}



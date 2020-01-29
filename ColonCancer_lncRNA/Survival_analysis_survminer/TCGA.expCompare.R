TCGA.expCompare <- function(cancerType = "COAD",geneName = "ENSG00000242268",
                            print = TRUE){
  #### 1.Read TCGA dataset ####
  filePath <- paste0("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/",cancerType,"_UCSC_Toil_tpm_dataset.rds")
  TCGA.Toil.tpm.dataset <- readRDS(filePath)
  # extract expression data
  TCGA.RSEM.gene.tpm <- TCGA.Toil.tpm.dataset[[paste0(cancerType,".RSEM.gene.tpm")]]
  # Remove the number after "." in Ensembl ID
  genelist <- rownames(TCGA.RSEM.gene.tpm)
  genelist <- gsub("\\..*","",genelist)
  rownames(TCGA.RSEM.gene.tpm) <- genelist
  #### judge the gene exist in genelist ####
  if(!(geneName %in% genelist)){
    print(paste0(geneName," not exist in TCGA data"))
  }else{
  # Zero percent in the data
  geneTPM <- TCGA.RSEM.gene.tpm[geneName,]
  nonZeroPercent <- round(sum(geneTPM != 0 ) / length(geneTPM), digits = 4)
  #### 2.Separete TCGA data into normal and cancer types
  sampleName <- colnames(TCGA.RSEM.gene.tpm)
  samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
  #table(samplType)
  ## "01" was primary tumor, "11' was normal tissue
  tumorIndex <- samplType == "01"
  normalIndex <- samplType == "11"
  ## Expression data for normal and cancer
  tumor.exp <- TCGA.RSEM.gene.tpm[geneName, tumorIndex]
  normal.exp<- TCGA.RSEM.gene.tpm[geneName, normalIndex]
  # Number of tumor and normal samples
  NumSam.tumor <- length(tumor.exp)
  NumSam.normal <- length(normal.exp)
  #### 3.Build ggplot dataframe
  df <- data.frame(expressionLog2_0.01 = c(tumor.exp,normal.exp),
                   TissueTypes = c(rep("Tumor",length(tumor.exp)),rep("Normal",length(normal.exp))))
  # log2(tpm + 0.01) transform the expression data
  df$expressionLog2_0.01 <- log2(df$expressionLog2_0.01+0.01)
  p_value <- t.test(df$expressionLog2_0.01 ~ df$TissueTypes )
  #p_value$p.value
  ## 4.Arrange ggplot and print the output
  ## ggplot 
  library(ggplot2)
  p <- ggplot(df, aes(TissueTypes,expressionLog2_0.01))
  if(print){
    res <- p + geom_boxplot(fill = c("dodgerblue3","firebrick3")) +  labs(x= paste0("Normal.N = ",NumSam.normal,
                                                                             " Tumor.N = ",NumSam.tumor),
                                                                   y = "Log2(tpm + 0.01) expression",
                                                                   title =paste0(geneName," in ", cancerType),
                                                                   subtitle = paste0("p-value = ", p_value$p.value))
    print(res)
    
  }else{
    res <- p + geom_boxplot(fill = c("dodgerblue3","firebrick3")) +  
                labs(x= paste0("Normal.N = ",NumSam.normal," Tumor.N = ",NumSam.tumor),
                     y = "Log2(tpm + 0.01) expression",
                     title =paste0(geneName," in ", cancerType),
                     subtitle = paste0("p-value = ", p_value$p.value))
    ggsave(paste0(geneName,"_in_", cancerType,"_expression_compare_analysis.pdf"), res,
           width = 6,height = 5)
  }
}
  
}

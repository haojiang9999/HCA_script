### Pre-build function for TCGA survival analysis ###
TCGA.surv <- function(cancerType = "COAD",geneName = "ENSG00000242268",
                      print = TRUE){
  #### 1.Read TCGA dataset ####
  #print = TRUE
  #geneName = "ENSG00000232956"
  filePath <- paste0("/data8t_4/JH/MyJobs/Read_dataset/UCSC_Toil/",cancerType,"_UCSC_Toil_tpm_dataset.rds")
  TCGA.Toil.tpm.dataset <- readRDS(filePath)
  TCGA.RSEM.gene.tpm <- TCGA.Toil.tpm.dataset[[paste0(cancerType,".RSEM.gene.tpm")]]
  TCGA.pheno <- TCGA.Toil.tpm.dataset[[paste0(cancerType,".pheno")]]
  # Remove the number after "." in Ensembl ID
  genelist <- rownames(TCGA.RSEM.gene.tpm)
  genelist <- gsub("\\..*","",genelist)
  rownames(TCGA.RSEM.gene.tpm) <- genelist
  #### judge the gene exist in genelist ####
  if(!(geneName %in% genelist)){
    print(paste0(geneName," not exist in TCGA data"))
  }else{
  #### 2. Extract target gene expression data ####
  #geneName <- "ENSG00000242268"
  geneTPM <- TCGA.RSEM.gene.tpm[geneName,]
  nonZeroPercent <- round(sum(geneTPM != 0 ) / length(geneTPM), digits = 4)
  #### 3. order samples by gene expression####
  orderedSamples <- names(geneTPM[order(geneTPM)])
  # Calculate sample gene expresasion quantiles
  library(dplyr)
  sampleQuartile <- ntile(geneTPM[order(geneTPM)], 4) 
  firstQ <- (sampleQuartile > 1)*1
  secQ <- (sampleQuartile > 2)*1
  thirdQ <- (sampleQuartile > 3)*1
  nonZeroQ <- (geneTPM[order(geneTPM)] != 0 )*1
  #table(nonZeroQ)
  ## 4. Build table for survival analysis
  surTime <- TCGA.pheno[orderedSamples,c("OS","OS.time")]
  TCGA.sur.df <- data.frame(surTime, firstQ = firstQ,
                            secQ = secQ, thirdQ = thirdQ,
                            nonZeroQ = nonZeroQ)
  ## 5. Select tumor samples
  ## sample type ###
  sampleName <- rownames(TCGA.sur.df)
  samplType<- unlist(lapply(strsplit(sampleName,".",fixed=TRUE), "[[", 4))
  table(samplType)
  ## "01" was primary tumor, "11' was normal tissue
  tumorIndex <- samplType == "01"
  #normalIndex <- samplType == "11"
  TCGA.sur.df <- TCGA.sur.df[tumorIndex,]
  # How many samples used in survival analysis?
  sampleNumber <- length(TCGA.sur.df$firstQ)
  # export the TCGA.sur.df
  TCGA.sur.df <<- TCGA.sur.df
  #### 6. Survival analysis plot ####
  library(survminer)
  require("survival")
  fit.firstQ <- survfit(Surv(OS.time, OS) ~ firstQ, data = TCGA.sur.df)
  fit.secQ <- survfit(Surv(OS.time, OS) ~ secQ, data = TCGA.sur.df)
  fit.thirdQ <- survfit(Surv(OS.time, OS) ~ thirdQ, data = TCGA.sur.df)
  fit.nonZeroQ <- survfit(Surv(OS.time, OS) ~ nonZeroQ, data = TCGA.sur.df)
  
  #### 6. List of ggsurvplots ####
  require("survminer")
  splots <- list()
  splots[[1]] <- ggsurvplot(fit.firstQ, pval = T,pval.method = T,
                            title = paste0(geneName," in ",cancerType," survival"),
                            legend.title = paste0("nonZeroPercent =", nonZeroPercent,
                                                  " N = ", sampleNumber))
  splots[[2]] <- ggsurvplot(fit.secQ,, pval = T,
                            title = paste0(geneName," in ",cancerType," survival"),
                            legend.title = paste0("nonZeroPercent =", nonZeroPercent,
                                                  " N = ", sampleNumber))
  splots[[3]] <- ggsurvplot(fit.thirdQ,  pval = T,
                            title = paste0(geneName," in ",cancerType," survival"),
                            legend.title = paste0("nonZeroPercent =", nonZeroPercent,
                                                  " N = ", sampleNumber))
  splots[[4]] <- ggsurvplot(fit.nonZeroQ, pval = T,
                            title = paste0(geneName," in ",cancerType," survival"),
                            legend.title = paste0("nonZeroPercent =", nonZeroPercent,
                                                  " N = ", sampleNumber))
  ## 7. Arrange multiple ggsurvplots and print the output
  if(print){
    print(arrange_ggsurvplots(splots, print = print,
                                ncol = 2, nrow = 2))
  }else{
    res <- arrange_ggsurvplots(splots, print = FALSE)
    ggsave(paste0(geneName,"_in_", cancerType,"_survival_analysis.pdf"), res,
           width = 12,height = 5)
  }
  }
  }


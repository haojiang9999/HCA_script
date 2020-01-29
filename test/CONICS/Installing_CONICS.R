### Install CONICS ####
install.packages("beanplot")
install.packages("mixtools")
install.packages("pheatmap")
install.packages("zoo")
install.packages("squash")
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
BiocManager::install("scran")
#install.packages("scran")
########## Important step in the install Rhdf5lib #######
install.packages('http://www.bioconductor.org/packages/release/bioc/src/contrib/Rhdf5lib_1.6.1.tar.gz', repos=NULL, type='source')
#install.packages("devtools")
devtools::install_github("diazlab/CONICS/CONICSmat", dep = FALSE)
devtools::install_github("diazlab/CONICS/CONICSmat", dep = FALSE)




# 
download.file(url = "https://singlecell.broadinstitute.org/single_cell/data/public/SCP12/oligodendroglioma-intra-tumor-heterogeneity?filename=OG_processed_data_portal.txt",
              destfile = "OG_processed_data_portal.txt")

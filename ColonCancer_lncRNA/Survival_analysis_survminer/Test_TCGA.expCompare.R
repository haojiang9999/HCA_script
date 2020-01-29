#### Test_TCGA.expCompare.R
TCGA.expCompare(cancerType = "COAD",geneName = "ENSG00000225972",
                print = TRUE)
res <- TCGA.expCompare(cancerType = "COAD",geneName = "ENSG00000225972",
                print = TRUE)
system.time(TCGA.expCompare(cancerType = "COAD",geneName = "ENSG00000242268",
                         print = TRUE))
TCGA.expCompare(cancerType = "COAD",geneName = "ENSG00000225972",
                print = F)




ggsave()
ggsave("test.pdf", res,
        width = 6,height = 5
       )

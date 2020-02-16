#### 2.test.R
ACC <- readRDS("ACC_UCSC_Toil_expected_count_dataset.rds")
test<- ACC$ACC.RSEM.gene.expected_count_round
test[1:5,1:5]

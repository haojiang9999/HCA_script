# separate samples by patients
sample.names <- as.character(sample.list[,1])
# how many cells in patients
table(unlist(lapply(strsplit(sample.names,"_"), '[[', 3)))
# for CRC01
sample.CRC01 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC01"]
# for CRC01
sample.CRC02 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC02"]
sample.CRC04 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC04"]
sample.CRC09 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC09"]
sample.CRC10 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC10"]
sample.CRC11 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC11"]
sample.CRC12 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC12"]
sample.CRC13 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC13"]
sample.CRC14 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC14"]
sample.CRC15 <- sample.names[unlist(lapply(strsplit(sample.names,"_"), '[[', 3)) == "CRC15"]

# cut the sample ID names to CRC13_LN2_377 type
CRC01.id <- sub("^[^C]*", "", sample.CRC01)
cbind(sample.CRC01,CRC01.id)
CRC02.id <- sub("^[^C]*", "", sample.CRC02)
CRC04.id <- sub("^[^C]*", "", sample.CRC04)
CRC09.id <- sub("^[^C]*", "", sample.CRC09)
CRC10.id <- sub("^[^C]*", "", sample.CRC10)
CRC11.id <- sub("^[^C]*", "", sample.CRC11)
CRC12.id <- sub("^[^C]*", "", sample.CRC12)
CRC13.id <- sub("^[^C]*", "", sample.CRC13)
CRC14.id <- sub("^[^C]*", "", sample.CRC14)
CRC15.id <- sub("^[^C]*", "", sample.CRC15)
save(CRC02.id, CRC04.id, CRC09.id, CRC10.id,CRC11.id,
     CRC12.id, CRC13.id, CRC14.id, CRC15.id, CRC01.id,
     sample.CRC01, sample.CRC02, sample.CRC04, sample.CRC09,
     sample.CRC10, sample.CRC11, sample.CRC12, sample.CRC13,
     sample.CRC14, sample.CRC15,
     file = "sample.id")
load(file = "test.id")

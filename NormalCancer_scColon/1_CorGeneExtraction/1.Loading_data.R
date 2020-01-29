#### Load data 
log2.CV1500.dataset <- readRDS("/data8t_4/JH/MyJobs/NormalCancer_scColon/log2.CV1500.dataset.rds")
Log2.expList <- log2.CV1500.dataset$Log2.expList
log2.scReference.list.CV.1500 <- log2.CV1500.dataset$log2.scReference.list.CV.1500
log2.Tang.colon.cancer.FPKM.500 <- Log2.expList$log2.Tang.colon.cancer.FPKM.500
Ref.Tang.Adult.colon <- log2.scReference.list.CV.1500$Tang.Adult.colon

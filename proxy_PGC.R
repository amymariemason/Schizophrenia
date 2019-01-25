rall <- read.delim("~/CVD_Schizophrenia/PGC/rall.txt")
PGC_snps<-rall[rall$snpid %in% all_SNPs$rsid,]
# 6 snps missing:  rs1998013, rs2247056, rs11246602, rs139385870, rs36083386, rs11442819  
# find proxies?


###############################


# compare to proxies
#rs1998013
proxy_rs1998013 <- read.delim("~/CVD_Schizophrenia/proxy_rs1998013.txt")
proxy_rs1998013$missing<-(proxy_rs1998013$RS_Number %in% rall$snpid)
proxy_rs1998013<-proxy_rs1998013[proxy_rs1998013$missing==TRUE,]
summary(proxy_rs1998013$Dprime)
summary(proxy_rs1998013$R2)


#rs2247056

proxy_rs2247056 <- read.delim("~/CVD_Schizophrenia/proxy_rs2247056.txt")
proxy_rs2247056$missing<-(proxy_rs2247056$RS_Number %in% rall$snpid)
proxy_rs2247056<-proxy_rs2247056[proxy_rs2247056$missing==TRUE,]
summary(proxy_rs2247056$Dprime)
summary(proxy_rs2247056$R2)


#rs11246602

proxy_rs11246602 <- read.delim("~/CVD_Schizophrenia/proxy_rs11246602.txt")
proxy_rs11246602$missing<-(proxy_rs11246602$RS_Number %in% rall$snpid)
proxy_rs11246602<-proxy_rs11246602[proxy_rs11246602$missing==TRUE,]
summary(proxy_rs11246602$Dprime)
# none in PGC

# rs139385870

proxy_rs139385870 <- read.delim("~/CVD_Schizophrenia/proxy_rs139385870.txt")
proxy_rs139385870$missing<-(proxy_rs139385870$RS_Number %in% rall$snpid)
proxy_rs139385870<-proxy_rs139385870[proxy_rs139385870$missing==TRUE,]
summary(proxy_rs139385870$Dprime)

# rs36083386

proxy_rs36083386 <- read.delim("~/CVD_Schizophrenia/proxy_rs36083386.txt")
proxy_rs36083386$missing<-(proxy_rs36083386$RS_Number %in% rall$snpid)
proxy_rs36083386<-proxy_rs36083386[proxy_rs36083386$missing==TRUE,]
summary(proxy_rs36083386$Dprime)


# rs11442819

proxy_rs11442819 <- read.delim("~/CVD_Schizophrenia/proxy_rs11442819.txt")
proxy_rs11442819$missing<-(proxy_rs11442819$RS_Number %in% rall$snpid)
proxy_rs11442819<-proxy_rs11442819[proxy_rs11442819$missing==TRUE,]
summary(proxy_rs11442819$Dprime)

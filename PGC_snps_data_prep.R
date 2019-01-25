####################################################################
# repeat with PGC
# author: amy mason
# date: jan 2019
##############################################################

# PGC schizophrenia outcomes

PGC_snps<-read.csv("~/CVD_Schizophrenia/inputs/PGC_cutdown.csv", header=TRUE)
names(PGC_snps)<-c("X", "chr", "rsid","Effect.allele_sch", "Association.alleles_sch", "pos",  "info", "or", "se_sch", "p_sch", "ngt")
PGC_snps$beta_sch<-log(PGC_snps$or)
keepnames<- c( "rsid","chr", "pos", "Effect.allele_sch", "Association.alleles_sch", "pos", "beta_sch","se_sch", "p_sch") 
PGC<-PGC_snps[,colnames(PGC_snps)%in%keepnames]
PGC$proxyby="no"

# proxy set matchings

proxy_PGC <- read.csv("~/CVD_Schizophrenia/proxy_PGC.csv")
PGC_snps2<-read.csv("~/CVD_Schizophrenia/inputs/PGC_cutdown_proxy.csv", header=TRUE)
names(PGC_snps2)<-c("X", "chr", "rsid","Effect.allele_sch", "Association.alleles_sch", "pos",  "info", "or", "se_sch", "p_sch", "ngt")
PGC_snps2$beta_sch<-log(PGC_snps2$or)
keepnames<- c( "rsid","chr", "pos", "Effect.allele_sch", "Association.alleles_sch", "pos", "beta_sch","se_sch", "p_sch") 
PGC2<-PGC_snps2[,colnames(PGC_snps)%in%keepnames]
PGC2$proxyby<-PGC2$rsid

# replace with the original snp data to allow matching
PGC_temp<-merge(PGC2, proxy_PGC,by.x="rsid", by.y="Proxy.snp")
as.character(PGC_temp$Effect.allele_sch)==as.character(PGC_temp$proxy_effect)
as.character(PGC_temp$Association.alleles_sch)==as.character(PGC_temp$proxy_other)
# 
keepnames<-c( "snp","chr.y", "pos.y", "snp_effect", "snp_ref", "beta_sch","se_sch", "p_sch", "proxyby") 
PGC2<-PGC_temp[,colnames(PGC_temp)%in%keepnames]
names(PGC2)<-c( "se_sch","p_sch", "beta_sch", "proxyby", "rsid","chr", "pos", "Effect.allele_sch", "Association.alleles_sch")
PGC$chr<-gsub("chr", "", PGC$chr)

# make columns numeric or characteristic

PGC[,c("chr", "pos","beta_sch", "se_sch","p_sch")] <- as.data.frame(sapply(PGC[,c("chr", "pos","beta_sch", "se_sch","p_sch")], as.numeric))
PGC$rsid<-as.character(PGC$rsid)
PGC$Effect.allele_sch<-as.character(PGC$Effect.allele_sch)
PGC$Association.alleles_sch<-as.character(PGC$Association.alleles_sch)
PGC$proxyby<-as.character(PGC$proxyby)

PGC2[,c("chr", "pos","beta_sch", "se_sch","p_sch")] <- as.data.frame(sapply(PGC2[,c("chr", "pos","beta_sch", "se_sch","p_sch")], as.numeric))
PGC2$rsid<-as.character(PGC2$rsid)
PGC2$Effect.allele_sch<-as.character(PGC2$Effect.allele_sch)
PGC2$Association.alleles_sch<-as.character(PGC2$Association.alleles_sch)
PGC2$proxyby<-as.character(PGC2$proxyby)


Schiz<-rbind(PGC, PGC2, stringsAsFactors=FALSE)

# ready to merge
write.csv(Schiz, file="~/CVD_Schizophrenia/inputs/PGC_plus_proxy.csv")

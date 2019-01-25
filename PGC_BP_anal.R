# Auther = Amy Mason
# Date = Nov 2018
# Purpose = merge UKBB data with exposure data 
#################################

Schiz_1 <- read.csv("~/CVD_Schizophrenia/inputs/PGC_plus_proxy.csv",stringsAsFactors=FALSE)                                                        

bp <- read.csv("~/CVD_Schizophrenia/data/bp_associations.csv", stringsAsFactors=FALSE)
names(bp)[4]<-"Effect.allele"
names(bp)[1]<-"chr"
names(bp)[2]<-"pos"
keepnames_bp<- c("rsid", "chr_pos",  "Effect.allele",  "trait", "beta" , "se" , "chr", "pos") 
bp2<-bp[, colnames(bp)%in% keepnames_bp]

names(bp2)<-c("chr", "pos",  "hg19_coordinates",  "Effect.allele","rsid",  "beta" , "se" , "trait") 



############## bp data merge
##### merge data on chromosome position

all_bp<-merge(bp2, Schiz_1, by=c("rsid"), all.x=TRUE, all.y=FALSE)
length(unique(all_bp[,"rsid"]))==nrow(bp)
# TRUE so merge is 1:1
all_bp$chr.x==all_bp$chr.y
all_bp$pos.x==all_bp$pos.y
# all true


### check alignment

all_bp$align<- ifelse((as.character(all_bp$Effect.allele)==as.character(all_bp$Effect.allele_sch))&(as.character(all_bp$Association.allele)==as.character(all_bp$Association.alleles_sch)), T, F)
summary(all_bp$align)
all_bp[all_bp$align==F,]
# all the problems are to do with insertions and deletions - rs36083386 are misaligned

## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- all_bp$rsid %in% c("rs36083386")
# do flip on schz results
all_bp$beta_sch[when.to.flip] <- - all_bp$beta_sch[when.to.flip]
New_effect<- all_bp[when.to.flip,"Effect.allele_sch"]
New_assos<- all_bp[when.to.flip,"Association.alleles_sch"]
all_bp[when.to.flip,"Effect.allele_sch"]<-New_assos
all_bp[when.to.flip,"Association.alleles_sch"]<-New_effect


#alignments now all match! Check for palindrome
all_bp$align3<- ifelse((all_bp$Effect.allele%in%c("A","T"))&(all_bp$Association.alleles_sch%in%c("A","T")), T, F)
all_bp$align3<- ifelse((all_bp$Effect.allele%in%c("C","G"))&(all_bp$Association.alleles_sch%in%c("C","G")), T, all_bp$align3)
summary(all_bp$align3)
bp_check<-all_bp[all_bp$align3==T,c("rsid", "beta", "beta_sch", "Effect.allele", "Effect.allele_sch", "Association.alleles_sch", "trait")]
bp_check[bp_check$trait!="PP",]
# checked beta values on phenoscanner for these 5: - agreement apart from rs2760061, rs9372498  & flipped for BP rs72812846 & both rs743757 [note only check SBP, DBP]
# do flip on schz results
when.to.flip1 <- all_bp$rsid %in% c("rs2760061", "rs9372498", "rs743757")
when.to.flip2<-all_bp$rsid %in% c( "rs72812846", "rs743757")
all_bp$beta_sch[when.to.flip] <- - all_bp$beta_sch[when.to.flip]
all_bp$beta[when.to.flip2] <- - all_bp$beta[when.to.flip2]



#### MR-model - DBP


mr.fit <- mr_input(DBP$beta, DBP$se, DBP$beta_sch, DBP$se_sch, snps=DBP$rsid)
mr_allmethods(mr.fit, method = "main")
p<-mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="bp & Schizophrenia")
p+xlab("DBP unit increase") + ylab("Schizophrenia PGC log(OR)")


#### MR-model - SBP


mr.fit <- mr_input(SBP$beta, SBP$se, SBP$beta_sch, SBP$se_sch, snps=SBP$rsid)
mr_allmethods(mr.fit, method = "main")
p<-mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="SBP & Schizophrenia")
p+xlab("SBP unit increase") + ylab("Schizophrenia PGC log(OR)")


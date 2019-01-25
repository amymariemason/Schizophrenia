# Auther = Amy Mason
# Date = Nov 2018
# Purpose = merge UKBB data with exposure data 
#################################

Schiz_1 <- read.csv("~/CVD_Schizophrenia/inputs/PGC_plus_proxy.csv",stringsAsFactors=FALSE)                                                        

#phenoscanner data for snp list for bmi
bmi2<-read.delim("~/CVD_Schizophrenia/data/bmi_rs_PhenoScanner.tsv",stringsAsFactors=FALSE)
names(bmi2)[5]<-"Effect.allele"
names(bmi2)[6]<-"Association.alleles"
bmi2<-bmi2[bmi2$pmid=="25673413",]
bmi2<-bmi2[bmi2$ancestry=="European",]
bmi2<-bmi2[bmi2$trait=="Body mass index",]
## this is now the same set as Steve's so use this. 


bmi2 <- cbind(bmi2, do.call("rbind", strsplit(as.character(bmi2$hg19_coordinates), ":")))
names(bmi2)[23]<-"chr"
names(bmi2)[24]<-"pos"
bmi2$chr<-gsub("chr", "", bmi2$chr)
keepnames_bmi<- c("rsid", "hg19_coordinates",  "Effect.allele",  "Association.alleles", "trait", "study"
                  , "pmid", "ancestry", "beta" , "se" , "p" ,  "direction" ,"n", "chr", "pos") 
bmi<-bmi2[, colnames(bmi2)%in% keepnames_bmi]

names(bmi)<-c("rsid", "hg19_coordinates",  "Effect.allele",  "Association.allele", "trait", "study",
              "pmid", "ancestry", "beta" , "se" , "p" ,  "direction" ,"n", "chr", "pos") 

############## bmi data merge
##### merge data on chromosome position

all_bmi<-merge(bmi, Schiz_1, by=c("rsid"), all.x=TRUE, all.y=FALSE)
length(unique(all_bmi[,"rsid"]))==nrow(bmi)
# TRUE so merge is 1:1
all_bmi$chr.x==all_bmi$chr.y
all_bmi$pos.x==all_bmi$pos.y
# all true


### check alignment

all_bmi$align<- ifelse((as.character(all_bmi$Effect.allele)==as.character(all_bmi$Effect.allele_sch))&(as.character(all_bmi$Association.allele)==as.character(all_bmi$Association.alleles_sch)), T, F)
summary(all_bmi$align)

## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- as.character(all_bmi$Effect.allele) != as.character(all_bmi$Effect.allele_sch)
# do flip on schz results
all_bmi$beta_sch[when.to.flip] <- - all_bmi$beta_sch[when.to.flip]
New_effect<- all_bmi[when.to.flip,"Effect.allele_sch"]
New_assos<- all_bmi[when.to.flip,"Association.alleles_sch"]
all_bmi[when.to.flip,"Effect.allele_sch"]<-New_assos
all_bmi[when.to.flip,"Association.alleles_sch"]<-New_effect

# check again
all_bmi$align2<- ifelse((all_bmi$Effect.allele==all_bmi$Effect.allele_sch)&(all_bmi$Association.allele==all_bmi$Association.alleles_sch), T, F)
summary(all_bmi$align2)

#alignments now all match! Check for palindrome
all_bmi$align3<- ifelse((all_bmi$Effect.allele%in%c("A","T"))&(all_bmi$Association.allele%in%c("A","T")), T, F)
all_bmi$align3<- ifelse((all_bmi$Effect.allele%in%c("C","G"))&(all_bmi$Association.allele%in%c("C","G")), T, all_bmi$align3)
summary(all_bmi$align3)
bmi_check<-all_bmi[all_bmi$align3==T,c("rsid", "beta", "beta_sch", "Effect.allele", "Association.allele", "Effect.allele_sch", "Association.alleles_sch", "trait")]
bmi_check
# checked beta values on phenoscanner for these 5: - agreement in all cases
# do flip on schz results


#### MR-model


mr.fit <- mr_input(all_bmi$beta, all_bmi$se, all_bmi$beta_sch, all_bmi$se_sch, snps=all_bmi$rsid)

mr_ivw(mr.fit)
mr_allmethods(mr.fit, method = "main")
p<-mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="BMI & Schizophrenia")
p+xlab("Body Mass Index IVNT") + ylab("Schizophrenia PGC log(OR)")

mr_plot(mr.fit, line="egger", orientate=T, interactive = F, labels=T)





## remove outlier and try again


mr.fit_o <- mr_input(all_bmi[all_bmi$rsid!="rs17024393",]$beta, all_bmi[all_bmi$rsid!="rs17024393",]$se, all_bmi[all_bmi$rsid!="rs17024393",]$beta_sch, all_bmi[all_bmi$rsid!="rs17024393",]$se_sch)
mr_allmethods(mr.fit_o, method = "main")
mr_plot(mr_allmethods(mr.fit_o, method = "main"))
mr_plot(mr.fit_o, line="egger", orientate=T, interactive = F, labels=T)

### clearly no relations between BMI and schitzophrenia

### leave-one-out

# lets play with mr_base

# exposure set

exposure<-all_bmi[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"bmi"
exp_dat <- format_data(exposure, type="exposure")

outcome<-all_bmi[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
names(outcome)<-c("SNP", "beta", "se", "effect_allele")
outcome$Phenotype<-"schizophrenia"
out_dat <- format_data(outcome, type="outcome")

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = out_dat,
  action =1 # as I have already manually checked the MAFs of palindromic SNPs manually
)

mr_heterogeneity(dat)
# no evidence of hetrogeneity

mr_pleiotropy_test(dat)
# no evidence of pleiotropy driving result

# forest plot of single snps
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

#forest plot of leave one out
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

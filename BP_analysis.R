# Auther = Amy Mason
# Date = Nov 2018
# Purpose = merge UKBB data with exposure data 
#################################

Schiz_sch_1 <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_1.out", header=TRUE, quote="\"")
keepnames<- c( "rsid","chromosome", "position", "alleleA", "alleleB",
              "all_maf", "frequentist_add_pvalue", "all_total", "frequentist_add_beta_1", "frequentist_add_se_1" )
# allele B is the effect allele
Schiz_1<-Schiz_sch_1[,colnames(Schiz_sch_1)%in%keepnames]
names(Schiz_1)<- c("rsid","chr", "pos", "Association.alleles_sch", "Effect.allele_sch",  
                      "n_sch",  "maf_sch" ,"p_sch" , "beta_sch", "se_sch" )                                                             

#steve data for snp list for bp

bp <- read.csv("~/CVD_Schizophrenia/data/bp_associations.csv")
names(bp)[4]<-"Effect.allele"
names(bp)[1]<-"chr"
names(bp)[2]<-"pos"
keepnames_bp<- c("rsid", "chr_pos",  "Effect.allele",  "trait", "beta" , "se" , "chr", "pos") 
bp2<-bp[, colnames(bp)%in% keepnames_bp]

names(bp2)<-c("chr", "pos",  "hg19_coordinates",  "Effect.allele","rsid",  "beta" , "se" , "trait") 

############## bmi data merge
##### merge data on chromosome position
missing<-as.data.frame("na")
all_bp<-merge(bp2, Schiz_1, by=c("rsid"), all.x=TRUE, all.y=FALSE)
length(unique(all_bp[,"rsid"]))==nrow(bp)
nrow(all_bp[all_bp$chr.x!=all_bp$chr.y,])
all_bp[is.na(all_bp$chr.y),c("rsid","trait", "chr.x","pos.x")]
# 9 variants missing from UK BB data
# discussed this with Steve - these 9 varients were also left out of the depression paper, so fine to ignore for this analysis
all_bp<-all_bp[!is.na(all_bp$chr.y),]
# 93 variants remaining
missing <- all_bp[is.na(all_bp$chr.y),c("rsid","trait", "chr.x","pos.x")]



### check alignment

all_bp$align<- ifelse((as.character(all_bp$Effect.allele)==as.character(all_bp$Effect.allele_sch)), T, F)
all_bp$align2<- ifelse((as.character(all_bp$Effect.allele)==as.character(all_bp$Association.alleles_sch)), T, F)
table(all_bp$align, all_bp$align2)

## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- as.character(all_bp$Effect.allele) == as.character(all_bp$Association.alleles_sch)
# do flip on schz results
all_bp$beta_sch[when.to.flip] <- - all_bp$beta_sch[when.to.flip]
New_effect<- all_bp[when.to.flip,"Effect.allele_sch"]
New_assos<- all_bp[when.to.flip,"Association.alleles_sch"]
all_bp[when.to.flip,"Effect.allele_sch"]<-New_assos
all_bp[when.to.flip,"Association.alleles_sch"]<-New_effect
  
# check again
all_bp$align2<- ifelse((as.character(all_bp$Effect.allele)==as.character(all_bp$Effect.allele_sch)), T, F)
summary(all_bp$align2)

#alignments now all match! Check for palindrome
all_bp$align3<- ifelse((all_bp$Effect.allele_sch%in%c("A","T"))&(all_bp$Association.alleles_sch%in%c("A","T")), T, F)
all_bp$align3<- ifelse((all_bp$Effect.allele_sch%in%c("C","G"))&(all_bp$Association.alleles_sch%in%c("C","G")), T, all_bp$align3)
summary(all_bp$align3)
bp_check<-all_bp[all_bp$align3==T,c("rsid", "maf_sch", "beta", "Effect.allele", "Association.alleles_sch", "trait")]
bp_check
# checked maf on phenoscanner for these 5: only rs4256980 is incorrectly aligned due to incorrect strand/palindrome problem
# do flip on schz results
when.to.flip <- all_bp$rsid%in%c("rs13112725", "rs2289081","rs2760061","rs72812846")
all_bp$beta_sch[when.to.flip] <- - all_bp$beta_sch[when.to.flip]


# seperate into sets by trait

DBP<-all_bp[all_bp$trait=="DBP",]
SBP<-all_bp[all_bp$trait=="SBP",]
PP<-all_bp[all_bp$trait=="PP",]


#### MR-model for DBP


mr.fit <- mr_input(DBP$beta, DBP$se, DBP$beta_sch, DBP$se_sch, snps=DBP$rsid)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main"))

# no significant fit


### clearly no relations between bp and schitzophrenia

### leave-one-out

# lets play with mr_base

# exposure set

exposure<-DBP[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"DBP"
exp_dat <- format_data(exposure, type="exposure")

outcome<-DBP[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
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


#### MR-model for SBP


mr.fit <- mr_input(SBP$beta, SBP$se, SBP$beta_sch, SBP$se_sch, snps=SBP$rsid)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main"))

# no significant fit


### clearly no relations between bp and schitzophrenia

### leave-one-out

# lets play with mr_base

# exposure set

exposure<-SBP[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"SBP"
exp_dat <- format_data(exposure, type="exposure")

outcome<-SBP[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
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
res_loo <- mr_leaveoneout(dat, method=mr_ivw)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]



#### MR-model for PP


mr.fit <- mr_input(PP$beta, PP$se, PP$beta_sch, PP$se_sch, snps=PP$rsid)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main"))

# no significant fit

# removing outliers


mr_plot(mr.fit, label=T, line="egger", orientate = T)
mr.fit <- mr_input(PP[PP$rsid!="rs79089478",]$beta, PP[PP$rsid!="rs79089478",]$se, PP[PP$rsid!="rs79089478",]$beta_sch, PP[PP$rsid!="rs79089478",]$se_sch, snps=PP[PP$rsid!="rs79089478",]$rsid)
mr_plot(mr.fit, label=T, line="egger", orientate = T)
mr_allmethods(mr.fit, method="main")


# exposure set

exposure<-PP[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"PP"
exp_dat <- format_data(exposure, type="exposure")

outcome<-PP[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
names(outcome)<-c("SNP", "beta", "se", "effect_allele")
outcome$Phenotype<-"schizophrenia"
out_dat <- format_data(outcome, type="outcome")

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = out_dat,
  action =1 # as I have already manually checked the MAFs of palindromic SNPs manually
)

mr_heterogeneity(dat)
# evidence of hetrogeneity!

mr_pleiotropy_test(dat)
# no evidence of pleiotropy driving result

# forest plot of single snps
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

#forest plot of leave one out
res_loo <- mr_leaveoneout(dat, method=mr_simple_median)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

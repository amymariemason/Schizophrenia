# Auther = Amy Mason
# Date = Nov 2018
# Purpose = merge UKBB data with exposure data 
#################################

library(MendelianRandomization)
library(epiDisplay)
library(TwoSampleMR)
library(ggplot2)


Schiz_sch_1 <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_1.out", header=TRUE, quote="\"")
keepnames<- c( "rsid","chromosome", "position", "alleleA", "alleleB",
               "all_maf", "frequentist_add_pvalue", "all_total", "frequentist_add_beta_1", "frequentist_add_se_1" )
# allele B is the effect allele
Schiz_1<-Schiz_sch_1[,colnames(Schiz_sch_1)%in%keepnames]
names(Schiz_1)<- c("rsid","chr", "pos", "Association.alleles_sch", "Effect.allele_sch",  
                   "n_sch",  "maf_sch" ,"p_sch" , "beta_sch", "se_sch" )                                                             

#steve data for snp list for bp

lipids<- read.csv("~/CVD_Schizophrenia/data/lipid_mrcat_all_MR_Catalogue_v2.csv")
lipids<-lipids[lipids$Source=="GLGC",]
lipids<-lipids[lipids$Year=="2010",]
names(lipids)[2]<-"rsid"
names(lipids)[3]<-"pos_name"

lipids <- cbind(lipids, do.call("rbind", strsplit(as.character(lipids$pos_name), ":")))
names(lipids)[25]<-"chr"
names(lipids)[26]<-"pos"
lipids$chr<-gsub("chr", "", lipids$chr)

keepnames_lipids<- c("rsid", "pos_name",  "Effect.Allele","Association.Alleles",  "Trait", "Beta" , "SE" , "chr", "pos", "MAF") 
lipids<-lipids[, colnames(lipids)%in% keepnames_lipids]

names(lipids)<-c("rsid", "hg19_coordinates","trait", "Effect.allele", "Association.alleles", "maf", "beta","se", "chr", "pos") 

############## data merge
##### merge data on chromosome position

all_lipids<-merge(lipids, Schiz_1, by=c("rsid"), all.x=TRUE, all.y=TRUE)
length(unique(all_lipids[,"rsid"]))==nrow(lipids)
nrow(all_lipids[all_lipids$chr.x!=all_lipids$chr.y,])
all_lipids[is.na(all_lipids$chr.y),]$hg19_coordinates
# 0 variants missing from UK BB data
all_lipids<-all_lipids[!is.na(all_lipids$chr.x),]



### check alignment

all_lipids$align<- ifelse((as.character(all_lipids$Effect.allele)==as.character(all_lipids$Effect.allele_sch)), T, F)
all_lipids$align2<- ifelse((as.character(all_lipids$Effect.allele)==as.character(all_lipids$Association.alleles_sch)), T, F)
# this table shows how many allele's don't match up 
table(all_lipids$align, all_lipids$align2)

# 1 variant (rs2954022) is not possible to align. lipid data gives A/C and Biobank only has C/T 
# - no way to change strand to match up
# DROP 1 var
all_lipids<-all_lipids[all_lipids$rsid!="rs2954022",]


## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- as.character(all_lipids$Effect.allele) == as.character(all_lipids$Association.alleles_sch)
# do flip on schz results
all_lipids$beta_sch[when.to.flip] <- - all_lipids$beta_sch[when.to.flip]
New_effect<- all_lipids[when.to.flip,"Effect.allele_sch"]
New_assos<- all_lipids[when.to.flip,"Association.alleles_sch"]
all_lipids[when.to.flip,"Effect.allele_sch"]<-New_assos
all_lipids[when.to.flip,"Association.alleles_sch"]<-New_effect

# check again
all_lipids$align2<- ifelse((as.character(all_lipids$Effect.allele)==as.character(all_lipids$Effect.allele_sch)), T, F)
summary(all_lipids$align2)

#alignments now all match! Check for palindrome
all_lipids$align3<- ifelse((all_lipids$Effect.allele_sch%in%c("A","T"))&(all_lipids$Association.alleles_sch%in%c("A","T")), T, F)
all_lipids$align3<- ifelse((all_lipids$Effect.allele_sch%in%c("C","G"))&(all_lipids$Association.alleles_sch%in%c("C","G")), T, all_lipids$align3)
summary(all_lipids$align3)
lipids_check<-all_lipids[all_lipids$align3==T,c("rsid", "maf_sch", "beta", "Effect.allele", "Effect.allele_sch", "trait", "maf_sch")]
lipids_check
# checked maf on phenoscanner for these 5: only rs4256980 is incorrectly aligned due to incorrect strand/palindrome problem
# no flipping needed



###################################################################### MR analysis

# seperate into sets by trait

HDL<-all_lipids[all_lipids$trait=="HDL",]
LDL<-all_lipids[all_lipids$trait=="LDL",]
Trigly<-all_lipids[all_lipids$trait=="Triglycerides",]

######################################################################

####HDL


#### MR-model for DBP


mr.fit <- mr_input(HDL$beta, HDL$se, HDL$beta_sch, HDL$se_sch, snps=HDL$rsid)
mr_allmethods(mr.fit, method = "main")
plot<-mr_plot(mr_allmethods(mr.fit, method = "main")) 
plot+ labs(title="HDL")

# no significant fit; check for outlier
mr_plot(mr.fit, labels=TRUE, interactive = FALSE) 

#HDL_out<-HDL[HDL$rsid!="",]
#mr.fit2 <- mr_input(HDL_out$beta, HDL_out$se, HDL_out$beta_sch, HDL_out$se_sch, snps=HDL_out$rsid)
#mr_plot(mr_allmethods(mr.fit2, method = "main")) 

### clearly no relations between bp and schitzophrenia

### leave-one-out

# lets play with mr_base

# exposure set

exposure<-HDL[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"HDL"
exp_dat <- format_data(exposure, type="exposure")

outcome<-HDL[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
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


##################################################################################
# LDL


mr.fit <- mr_input(LDL$beta, LDL$se, LDL$beta_sch, LDL$se_sch, snps=LDL$rsid)
mr_allmethods(mr.fit, method = "main")
plot<-mr_plot(mr_allmethods(mr.fit, method = "main")) 
plot+ labs(title="LDL")

# no significant fit; check for outlier
mr_plot(mr.fit, labels=TRUE, interactive = FALSE, line="IVW") 

LDL_out<-LDL[LDL$rsid!="rs646776",]
mr.fit2 <- mr_input(LDL_out$beta, LDL_out$se, LDL_out$beta_sch, LDL_out$se_sch, snps=LDL_out$rsid)
mr_plot(mr_allmethods(mr.fit2, method = "ivw")) 
mr_plot(mr.fit2, labels=TRUE, interactive = FALSE, line="ivw") +labs(title="LDL; re6776 excluded; ivw")

### clearly no relations between bp and schitzophrenia

### leave-one-out

# lets play with mr_base

# exposure set

exposure<-LDL_out[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"LDL"
exp_dat <- format_data(exposure, type="exposure")

outcome<-LDL_out[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
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
p3[[1]]+labs(title="LDL; re6776 excluded")

##################################################################################
# trigylerides


mr.fit <- mr_input(Trigly$beta, Trigly$se, Trigly$beta_sch, Trigly$se_sch, snps=Trigly$rsid)
mr_allmethods(mr.fit, method = "main")
plot<-mr_plot(mr_allmethods(mr.fit, method = "main")) 
plot+ labs(title="Triglyerides")

# no significant fit; check for outlier
mr_plot(mr.fit, labels=TRUE, interactive = FALSE) 

Trigly_out<-Trigly[Trigly$rsid!="rs1998013",]
mr.fit2 <- mr_input(Trigly_out$beta, Trigly_out$se, Trigly_out$beta_sch, Trigly_out$se_sch, snps=Trigly_out$rsid)
mr_plot(mr_allmethods(mr.fit2, method = "ivw")) 
mr_plot(mr.fit2, labels=TRUE, interactive = FALSE, line="ivw") +labs(title="Trigly; re6776 excluded; ivw")

### clearly no relations between bp and schitzophrenia

### leave-one-out

# lets play with mr_base

# exposure set

exposure<-Trigly[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"Trigly"
exp_dat <- format_data(exposure, type="exposure")

outcome<-Trigly[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
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


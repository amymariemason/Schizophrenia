# Auther = Amy Mason
# Date = Jan 2019
# Purpose = merge UKBB data with exposure data 
#################################

library(MendelianRandomization)
library(epiDisplay)
library(TwoSampleMR)
library(ggplot2)
library(readxl)

#######
#DATA MATCH
#######

Schiz_sch_1 <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_2.out", header=TRUE, quote="\"")
keepnames<- c( "rsid","chromosome", "position", "alleleA", "alleleB",
               "all_maf", "frequentist_add_pvalue", "all_total", "frequentist_add_beta_1", "frequentist_add_se_1" )
# allele B is the effect allele
Schiz_1<-Schiz_sch_1[,colnames(Schiz_sch_1)%in%keepnames]
names(Schiz_1)<- c("rsid","chr", "pos", "Association.alleles_sch", "Effect.allele_sch",  
                   "n_sch",  "maf_sch" ,"p_sch" , "beta_sch", "se_sch" )                                                             


# Inflammatory snps
inflamatory<- read_excel("~/CVD_Schizophrenia/data/inflamatory_markers_v2.xlsx")
inflamatory[7,1]<-"CRP"
inflamatory[13:15,2]<-"ICAM-1"
inflamatory[6,2]<-"TNF"

inflamatory <- cbind(inflamatory, do.call("rbind", strsplit(as.character(inflamatory$`chr:locus`), ":")))
names(inflamatory)[11]<-"chr"
names(inflamatory)[12]<-"pos"
inflamatory$chr<-gsub("chr", "", inflamatory$chr)

# match names to extract
keepnames<- c("locus", "inflammatory.factors",  "effect.allele","reference allele",  "beta" , "se" , "chr", "pos", "MAF") 
inflamatory<-inflamatory[, colnames(inflamatory)%in% keepnames]
names(inflamatory)<-c("trait", "rsid", "Effect.allele", "Association.alleles", "beta", "maf","se", "chr", "pos") 


##### merge data on chromosome position

all_inflamatory<-merge(inflamatory, Schiz_1, by=c("rsid"), all.x=TRUE, all.y=FALSE)
length(unique(all_inflamatory[,"rsid"]))==nrow(inflamatory)
nrow(inflamatory[inflamatory$chr.x!=inflamatory$chr.y,])
all_inflamatory[is.na(all_inflamatory$chr.y),c("rsid","trait", "chr.x","pos.x")]
# 1 variants missing from UK BB data - this is resolved with additional dataset at end
all_inflamatory<-all_inflamatory[!is.na(all_inflamatory$chr.y),]

#################################################################################################
### check alignment
################################################################################################

all_inflamatory$align<- ifelse((as.character(all_inflamatory$Effect.allele)==as.character(all_inflamatory$Effect.allele_sch)), T, F)
all_inflamatory$align2<- ifelse((as.character(all_inflamatory$Effect.allele)==as.character(all_inflamatory$Association.alleles_sch)), T, F)
# this table shows how many allele's don't match up 
table(all_inflamatory$align, all_inflamatory$align2)
# 3 with definite strand issues, investigate
# label as if different strand
when.to.flip <- all_inflamatory$align == F & all_inflamatory$align2==F
all_inflamatory$effect_temp<-all_inflamatory$Effect.allele_sch

#flip alleles clearly on wrong strand
all_inflamatory[when.to.flip & all_inflamatory$Effect.allele_sch=="T",]$effect_temp<-"A"
#all_inflamatory[when.to.flip & all_inflamatory$Effect.allele_sch=="A",]$effect_temp<-"T"
#all_inflamatory[when.to.flip & all_inflamatory$Effect.allele_sch=="C",]$effect_temp<-"G"
all_inflamatory[when.to.flip & all_inflamatory$Effect.allele_sch=="G",]$effect_temp<-"C"
#
all_inflamatory$ass_temp<-all_inflamatory$Association.alleles_sch
all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="T",]$ass_temp<-"A"
#all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="A",]$ass_temp<-"T"
all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="C",]$ass_temp<-"G"
#all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="G",]$ass_temp<-"C"

all_inflamatory$Effect.allele_sch<-all_inflamatory$effect_temp
all_inflamatory$Association.alleles_sch<-all_inflamatory$ass_temp

# check again

all_inflamatory$align<- ifelse((as.character(all_inflamatory$Effect.allele)==as.character(all_inflamatory$Effect.allele_sch)), T, F)
all_inflamatory$align2<- ifelse((as.character(all_inflamatory$Effect.allele)==as.character(all_inflamatory$Association.alleles_sch)), T, F)
table(all_inflamatory$align, all_inflamatory$align2)



# now all alleles match on strand, unless palindrome

# sort mis-alignments

## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- as.character(all_inflamatory$Effect.allele) == as.character(all_inflamatory$Association.alleles_sch)
# do flip on schz results
all_inflamatory$beta_sch[when.to.flip] <- - all_inflamatory$beta_sch[when.to.flip]
New_effect<- all_inflamatory[when.to.flip,"Effect.allele_sch"]
New_assos<- all_inflamatory[when.to.flip,"Association.alleles_sch"]
all_inflamatory[when.to.flip,"Effect.allele_sch"]<-New_assos
all_inflamatory[when.to.flip,"Association.alleles_sch"]<-New_effect

# check again
all_inflamatory$align2<- ifelse((as.character(all_inflamatory$Effect.allele)==as.character(all_inflamatory$Effect.allele_sch)), T, F)
summary(all_inflamatory$align2)


# check palindromes

all_inflamatory$align3<- ifelse((all_inflamatory$Effect.allele_sch%in%c("A","T"))&(all_inflamatory$Association.alleles_sch%in%c("A","T")), T, F)
all_inflamatory$align3<- ifelse((all_inflamatory$Effect.allele_sch%in%c("C","G"))&(all_inflamatory$Association.alleles_sch%in%c("C","G")), T, all_inflamatory$align3)
summary(all_inflamatory$align3)
inflamatory_check<-all_inflamatory[all_inflamatory$align3==T,c("rsid","maf", "maf_sch", "beta", "Effect.allele", "Association.alleles", "Effect.allele_sch", "Association.alleles_sch", "trait")]
inflamatory_check
# from maf phenoscanner results, clearly mislabelled on 1 variant - flip that
# do flip on schz results
when.to.flip <- all_inflamatory$rsid%in%c("rs1800947")
all_inflamatory$beta_sch[when.to.flip] <- - all_inflamatory$beta_sch[when.to.flip]


all_inflamatory$beta<-as.numeric(as.character(all_inflamatory$beta))
all_inflamatory$se<-as.numeric(as.character(all_inflamatory$se))

# seperate into sets by trait

CRP<-all_inflamatory[all_inflamatory$trait=="CRP",]
IL1<-all_inflamatory[all_inflamatory$trait=="IL1",]
IL6<-all_inflamatory[all_inflamatory$trait=="IL6",]
P_sel<-all_inflamatory[all_inflamatory$trait=="P-selectin",]
SICAM<-all_inflamatory[all_inflamatory$trait=="sICAM-1",]
TNF<-all_inflamatory[all_inflamatory$trait=="TNF-alpha",]


########################################################################################################################
# ANALYSIS - multiple snps
########################################################################################################################

###CRP


mr.fit <- mr_input(CRP$beta, CRP$se, CRP$beta_sch, CRP$se_sch, snps=CRP$rsid)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="CRP")

# MR-base 

# exposure set

exposure<-CRP[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"CRP"
exp_dat <- format_data(exposure, type="exposure")

outcome<-CRP[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
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


### IL1



mr.fit <- mr_input(IL1$beta, IL1$se, IL1$beta_sch, IL1$se_sch, snps=IL1$rsid)
mr_ivw(mr.fit)
mr_plot(mr.fit, interactive=F) +labs(title="IL1")




### IL6


mr.fit <- mr_input(IL6$beta, IL6$se, IL6$beta_sch, IL6$se_sch, snps=IL6$rsid)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="IL6")

## sICAM

mr.fit <- mr_input(SICAM$beta, SICAM$se, SICAM$beta_sch, SICAM$se_sch, snps=SICAM$rsid)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="SICAM")

############################################################################################################
# single snp analysis
############################################################################################################

sch <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_1.out", header=TRUE, quote="\"")
#* P-selectin
# increase per allele copy

snpid = "rs6136"

round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]), 2)
round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]-1.96*sch[which(sch$rsid==snpid), "frequentist_add_se_1"]), 2)
round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]+1.96*sch[which(sch$rsid==snpid), "frequentist_add_se_1"]), 2)
sch[which(sch$rsid==snpid), "frequentist_add_pvalue"]


#* Fibrinogen

sch <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_1_repeats.out", header=TRUE, quote="\"") 

snpid = "rs7439150"

round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]), 2)
round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]-1.96*sch[which(sch$rsid==snpid), "frequentist_add_se_1"]), 2)
round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]+1.96*sch[which(sch$rsid==snpid), "frequentist_add_se_1"]), 2)
sch[which(sch$rsid==snpid), "frequentist_add_pvalue"]


# TNF-alpha


snpid = "rs1800629"

round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]), 2)
round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]-1.96*sch[which(sch$rsid==snpid), "frequentist_add_se_1"]), 2)
round(exp(sch[which(sch$rsid==snpid), "frequentist_add_beta_1"]+1.96*sch[which(sch$rsid==snpid), "frequentist_add_se_1"]), 2)
sch[which(sch$rsid==snpid), "frequentist_add_pvalue"]

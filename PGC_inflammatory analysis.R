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


Schiz_1 <- read.csv("~/CVD_Schizophrenia/inputs/PGC_plus_proxy.csv",stringsAsFactors=FALSE)                                                         


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
# all there

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
#all_inflamatory[when.to.flip & all_inflamatory$Effect.allele_sch=="G",]$effect_temp<-"C"
#
all_inflamatory$ass_temp<-all_inflamatory$Association.alleles_sch
#all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="T",]$ass_temp<-"A"
#all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="A",]$ass_temp<-"T"
all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="C",]$ass_temp<-"G"
all_inflamatory[when.to.flip & all_inflamatory$Association.alleles_sch=="G",]$ass_temp<-"C"

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
inflamatory_check<-all_inflamatory[all_inflamatory$align3==T,c("rsid","maf", "beta", "Effect.allele", "Association.alleles", "Effect.allele_sch", "Association.alleles_sch", "trait")]
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
Fibrinogen<- all_inflamatory[all_inflamatory$trait=="Fibrinogen",]

########################################################################################################################
# ANALYSIS - multiple snps
########################################################################################################################

###CRP

library(TwoSampleMR)
rho<-ld_matrix(CRP$rsid)
detach("package:TwoSampleMR", unload=TRUE)

# directions of snps differs
rho[c(1,4),2:3]<- -rho[c(1,4),2:3]
rho[2:3, c(1,4)]<- -rho[2:3, c(1,4)]

mr.fit <- mr_input(CRP$beta, CRP$se, CRP$beta_sch, CRP$se_sch, snps=CRP$rsid, correlation=rho)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="CRP")
mr_plot(mr.fit, interactive=FALSE, method="mr_egger") +labs(title="CRP")

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

library(TwoSampleMR)
rho<-ld_matrix(IL1$rsid)
detach("package:TwoSampleMR", unload=TRUE)

# directions of snps differ
rho[2,1]<--rho[2,1]
rho[1,2]<--rho[1,2]

mr.fit <- mr_input(IL1$beta, IL1$se, IL1$beta_sch, IL1$se_sch, snps=IL1$rsid, correlation=rho)
mr_ivw(mr.fit)
mr_plot(mr.fit, interactive=F) +labs(title="IL1")




### IL6
library(TwoSampleMR)
rho<-ld_matrix(IL6$rsid)
detach("package:TwoSampleMR", unload=TRUE)

# directions of snps differ
rho[1:2,3]<--rho[1:2,3]
rho[3,1:2]<--rho[3, 1:2]


mr.fit <- mr_input(IL6$beta, IL6$se, IL6$beta_sch, IL6$se_sch, snps=IL6$rsid, correlation=rho)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="IL6")

## sICAM
library(TwoSampleMR)
rho<-ld_matrix(SICAM$rsid)
detach("package:TwoSampleMR", unload=TRUE)

# directions of snps differ
rho
SICAM[,c("rsid", "Effect.allele", "Association.alleles")]
#no changes needed

mr.fit <- mr_input(SICAM$beta, SICAM$se, SICAM$beta_sch, SICAM$se_sch, snps=SICAM$rsid,correlation=rho)
mr_allmethods(mr.fit, method = "main")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="SICAM")

############################################################################################################
# single snp analysis
############################################################################################################

sch <- read.csv("~/CVD_Schizophrenia/inputs/PGC_plus_proxy.csv",stringsAsFactors=FALSE)      
#* P-selectin
# increase per allele copy

snpid = "rs6136"

round(exp(sch[which(sch$rsid==snpid), "beta_sch"]), 2)
round(exp(sch[which(sch$rsid==snpid), "beta_sch"]-1.96*sch[which(sch$rsid==snpid), "se_sch"]), 2)
round(exp(sch[which(sch$rsid==snpid), "beta_sch"]+1.96*sch[which(sch$rsid==snpid), "se_sch"]), 2)
sch[which(sch$rsid==snpid), "beta_sch"]



exposure<-P_sel[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"P_selectin"
exp_dat <- format_data(exposure, type="exposure")

outcome<-P_sel[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
names(outcome)<-c("SNP", "beta", "se", "effect_allele")
outcome$Phenotype<-"schizophrenia"
out_dat <- format_data(outcome, type="outcome")

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = out_dat,
  action =1 # as I have already manually checked the MAFs of palindromic SNPs manually
)

#Calculate ration
keep<-mr(dat)
exp(keep$b)
exp(keep$b)+1.96*keep$se
 exp(keep$b)-1.96*keep$se

#* Fibrinogen

snpid = "rs7439150"

round(exp(sch[which(sch$rsid==snpid), "beta_sch"]), 2)
round(exp(sch[which(sch$rsid==snpid), "beta_sch"]-1.96*sch[which(sch$rsid==snpid), "se_sch"]), 2)
round(exp(sch[which(sch$rsid==snpid), "beta_sch"]+1.96*sch[which(sch$rsid==snpid), "se_sch"]), 2)
sch[which(sch$rsid==snpid), "beta_sch"]



exposure<-Fibrinogen[,c("rsid", "beta", "se","Effect.allele")]
names(exposure)<-c("SNP", "beta", "se", "effect_allele")
exposure$Phenotype<-"Fibrinogen"
exp_dat <- format_data(exposure, type="exposure")

outcome<-Fibrinogen[,c("rsid", "beta_sch", "se_sch","Effect.allele_sch")]
names(outcome)<-c("SNP", "beta", "se", "effect_allele")
outcome$Phenotype<-"schizophrenia"
out_dat <- format_data(outcome, type="outcome")

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = out_dat,
  action =1 # as I have already manually checked the MAFs of palindromic SNPs manually
)



# TNF-alpha


snpid = "rs1800629"

round(exp(sch[which(sch$rsid==snpid), "beta_sch"]), 2)
round(exp(sch[which(sch$rsid==snpid), "beta_sch"]-1.96*sch[which(sch$rsid==snpid), "se_sch"]), 2)
round(exp(sch[which(sch$rsid==snpid), "beta_sch"]+1.96*sch[which(sch$rsid==snpid), "se_sch"]), 2)
sch[which(sch$rsid==snpid), "beta_sch"]




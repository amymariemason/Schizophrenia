################################################################
# Larger fibrinogen panel compared with depression and schizophrenia
# author: amy mason
# date: jan 2019
# look at set of snps from phenoscanner; subset of those from FGA or FGB
################################################################
library(MendelianRandomization)
library(epiDisplay)
library(TwoSampleMR)
library(ggplot2)

##


Schiz_sch_1 <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_1_fibro.out", header=TRUE, quote="\"")
keepnames<- c( "rsid","chromosome", "position", "alleleA", "alleleB",
               "all_maf", "frequentist_add_pvalue", "all_total", "frequentist_add_beta_1", "frequentist_add_se_1" )
# allele B is the effect allele
Schiz_1<-Schiz_sch_1[,colnames(Schiz_sch_1)%in%keepnames]
names(Schiz_1)<- c("rsid","chr", "pos", "Association.alleles_sch", "Effect.allele_sch",  
                   "n_sch",  "maf_sch" ,"p_sch" , "beta_sch", "se_sch" )  

# include original snp result


sch <- read.table("~/CVD_Schizophrenia/inputs/Schiz_sch_1_repeats.out", header=TRUE, quote="\"") 
keepnames<- c( "rsid","chromosome", "position", "alleleA", "alleleB",
               "all_maf", "frequentist_add_pvalue", "all_total", "frequentist_add_beta_1", "frequentist_add_se_1" )
# allele B is the effect allele
Sch<-sch[,colnames(sch)%in%keepnames]
names(Sch)<- c("rsid","chr", "pos", "Association.alleles_sch", "Effect.allele_sch",  
                   "n_sch",  "maf_sch" ,"p_sch" , "beta_sch", "se_sch" )  

schiz<-rbind(Schiz_1,Sch)

##################################

# phenoscanner results

# Inflammatory snps (phenoscanner)
Fibrinogen_PhenoScanner_GWAS <- read.delim("~/CVD_Schizophrenia/data/Fibrinogen_PhenoScanner_GWAS.tsv")
Fibri <-Fibrinogen_PhenoScanner_GWAS[!is.na(Fibrinogen_PhenoScanner_GWAS$beta),] 
Fibri <-Fibri[Fibri$ancestry=="European",] 
Fibri$wanted <-as.character(Fibri$hgnc) %in% c("FGA", "FGB", "FGG")
Fibri<-Fibri[Fibri$wanted==T,]
Fibri <-Fibri[Fibri$p<5*10^(-8),] 
Fibri$unit<-as.character(Fibri$unit)
# fix unit
Fibri[Fibri$unit=="mg/dl increase",]$beta<-as.numeric(Fibri[Fibri$unit=="mg/dl increase",]$beta)/100
Fibri[Fibri$unit=="mg/dl increase",]$unit<-"g/L increase"

Fibri <- cbind(Fibri, do.call("rbind", strsplit(as.character(Fibri$hg19_coordinates), ":")))
names(Fibri)[33:34]<-c("chr","pos")
Fibri$chr<-gsub("chr","",Fibri$chr)
 
keepnames<- c("rsid", "a1", "a2", "trait", "beta", "se", "p", "chr", "pos")
Fibri<-Fibri[, colnames(Fibri)%in% keepnames]
names(Fibri)<- c("rsid", "Effect.allele", "Association.alleles", "trait", "beta", "se","p" , "chr", "pos") 
Fibri$maf<-"na"


# Inflammatory snps ()
inflamatory<- read_excel("~/CVD_Schizophrenia/data/inflamatory_markers_v2.xlsx")
inflamatory <- cbind(inflamatory, do.call("rbind", strsplit(as.character(inflamatory$`chr:locus`), ":")))
names(inflamatory)[11]<-"chr"
names(inflamatory)[12]<-"pos"
inflamatory$chr<-gsub("chr", "", inflamatory$chr)
inflamatory<-inflamatory[inflamatory$locus=="Fibrinogen",]

keepnames<- c("locus", "inflammatory.factors",  "effect.allele","reference allele",  "beta" , "se" , "chr", "pos", "MAF") 
inflamatory<-inflamatory[, colnames(inflamatory)%in% keepnames]
names(inflamatory)<-c("trait", "rsid", "Effect.allele", "Association.alleles", "beta", "maf","se", "chr", "pos")
inflamatory$p<-""

inflam<-rbind(Fibri, inflamatory)


##### merge data on chromosome position

all_Fibri<-merge(inflam, schiz, by=c("rsid"), all.x=TRUE, all.y=FALSE)
length(unique(all_Fibri[,"rsid"]))==nrow(inflam)
nrow(all_Fibri[all_Fibri$chr.x!=all_Fibri$chr.y,])
# full match

#################################################################################################
### check alignment
################################################################################################

all_Fibri$align<- ifelse((as.character(all_Fibri$Effect.allele)==as.character(all_Fibri$Effect.allele_sch)), T, F)
all_Fibri$align2<- ifelse((as.character(all_Fibri$Effect.allele)==as.character(all_Fibri$Association.alleles_sch)), T, F)
# this table shows how many allele's don't match up 
table(all_Fibri$align, all_Fibri$align2)
# 0 with obvious strand issues
# now all alleles match on strand, unless palindrome

# sort mis-alignments

## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- as.character(all_Fibri$Effect.allele) == as.character(all_Fibri$Association.alleles_sch)
# do flip on schz results
all_Fibri$beta_sch[when.to.flip] <- - all_Fibri$beta_sch[when.to.flip]
New_effect<- all_Fibri[when.to.flip,"Effect.allele_sch"]
New_assos<- all_Fibri[when.to.flip,"Association.alleles_sch"]
all_Fibri[when.to.flip,"Effect.allele_sch"]<-New_assos
all_Fibri[when.to.flip,"Association.alleles_sch"]<-New_effect

# check again
all_Fibri$align2<- ifelse((as.character(all_Fibri$Effect.allele)==as.character(all_Fibri$Effect.allele_sch)), T, F)
summary(all_Fibri$align2)


# check palindromes

all_Fibri$align3<- ifelse((all_Fibri$Effect.allele_sch%in%c("A","T"))&(all_Fibri$Association.alleles_sch%in%c("A","T")), T, F)
all_Fibri$align3<- ifelse((all_Fibri$Effect.allele_sch%in%c("C","G"))&(all_Fibri$Association.alleles_sch%in%c("C","G")), T, all_Fibri$align3)
summary(all_Fibri$align3)
inflamatory_check<-all_Fibri[all_Fibri$align3==T,c("rsid","maf", "maf_sch", "beta", "Effect.allele", "Association.alleles", "Effect.allele_sch", "Association.alleles_sch", "trait")]
inflamatory_check
# from maf phenoscanner results, checked and fine



all_Fibri$beta<-as.numeric(as.character(all_Fibri$beta))
all_Fibri$se<-as.numeric(as.character(all_Fibri$se))
all_Fibri$beta_sch<-as.numeric(as.character(all_Fibri$beta_sch))
all_Fibri$se_sch<-as.numeric(as.character(all_Fibri$se_sch))

# R2 - there is correlation between many of these snps

all_Fibri_LD<-all_Fibri[all_Fibri$rsid %in% c("rs7439150","rs1800788", "rs6050"),]
R2 <- read.delim("~/CVD_Schizophrenia/data/fibrinogen_R^2.txt")
R2<-as.matrix(R2[1:5,2:6])
rownames(R2)<-colnames(R2)

all_Fibri$reference <- sapply(all_Fibri$rsid, function(x) which(x == rownames(R2)))
keep_rows<-all_Fibri[all_Fibri$rsid %in% c("rs7439150","rs1800788", "rs6050"),]$reference
R2_LD<-R2[keep_rows,keep_rows]

###########################################################################
# ANALYSIS


mr.fit <- mr_input(all_Fibri_LD$beta, all_Fibri_LD$se, all_Fibri_LD$beta_sch, all_Fibri_LD$se_sch, snps=all_Fibri_LD$rsid, correlation = as.matrix(R2_LD))
mr_allmethods(mr.fit, method="ivw")
mr_plot(mr_allmethods(mr.fit, method = "main")) +labs(title="Fibrinogen & Schizophrenia")
+xlab="Fibrinogen g/L" + ylab="Schizophrenia UKBB log(OR)"


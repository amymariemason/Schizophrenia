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




# fasting glucose
fasting_glucose <- read.csv("~/CVD_Schizophrenia/data/fasting_glucose_assoc.csv")
fasting_insulin <- read.csv("~/CVD_Schizophrenia/data/fasting_insulin_assoc.csv")
HOMA_IR <- read.csv("~/CVD_Schizophrenia/data/HOMA_IR_assoc.csv")
leptin<- read.csv("~/CVD_Schizophrenia/data/leptin_assoc.csv")

# stick these together for ease of merging with schizophrenia data


all_gluc<-rbind(fasting_glucose,fasting_insulin)
names(all_gluc)[5]="Effect.allele"
names(all_gluc)[6]="Association.allele"

# add chr/pos
all_gluc <- cbind(all_gluc, do.call("rbind", strsplit(as.character(all_gluc$hg19_coordinates), ":")))
names(all_gluc)[16]<-"chr"
names(all_gluc)[17]<-"pos"
all_gluc$chr<-gsub("chr", "", all_gluc$chr)

# remove unneded columns and convert to character
all_gluc<-all_gluc[,c("rsid", "Effect.allele", "Association.allele", "trait",  "beta", "se","p","direction", "chr", "pos" )]
i <- sapply(all_gluc, is.factor)
all_gluc[i] <- lapply(all_gluc[i], as.character)



# add HOMA_IR
names(HOMA_IR)<-c("rsid", "Effect.allele", "Association.allele", "beta", "se","p", "direction", "chr", "pos")
HOMA_IR$trait<-"HOMA_IR"
i <- sapply(HOMA_IR, is.factor)
HOMA_IR[i] <- lapply(HOMA_IR[i], as.character)

all_gluc<-rbind(all_gluc,HOMA_IR)


# add leptin
leptin<-leptin[,c("rsid","chr", "pos","effect", "other", "beta", "se", "p")]
names(leptin)<-c("rsid","chr", "pos","Effect.allele", "Association.allele", "beta", "se", "p")
leptin$direction<-NA
leptin$trait<-"leptin"

all_gluc<-rbind(all_gluc,leptin)

########## merge with schizophrenia data
##### merge data on chromosome position

all<-merge(all_gluc, Schiz_1, by=c("rsid"), all.x=TRUE, all.y=FALSE)
length(unique(all[,"rsid"]))==nrow(all_gluc)
nrow(all[all$chr.x!=all$chr.y,])
all[is.na(all$chr.y),c("rsid","trait", "chr.x","pos.x")]
# 1 variants missing from UK BB data
#AMY NOTE: FIX MISSING VARIANT
all<-all[!is.na(all$chr.y),]




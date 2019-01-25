#
# Written by Amy Mason
# Part of investigation into effect of history of CVD 
# and other co-morbidity factors on schizophrenia
#
# This file chooses SNP sets for insulin and glucose related snp sets, 
# using data from referenced papers and from Phenoscanner results using the
# highlighed key word and a p-value cut-off of 10E-8. Where varients are reported
# multiple times, the value from the large set of cases/controls is used.
#
#
library(dplyr)
library(tidyr)

setwd("U://My Documents/CVD_Schizophrenia/")

# functions


qctoolify <- function(table, chr=chr, pos=pos, filename){
  #table is the table of chromosomes and positions you want to extract
  #chr is the name of the column containing chromosone information (integrer)
  #pos is the name of the column containing position information
  #filename is where the file should be saved
  ###################
  # this function creates a space seperated list of snps ready to use in qctool
  table$output<-ifelse(as.numeric(table$chr)<10, paste0("0",table$chr,":", table$pos), paste0(table$chr,":", table$pos))
  write.table(t(table$output),filename,sep=" ",row.names=FALSE, quote =FALSE, col.names=FALSE)
}


################################################################################################################################
# fasting glucose snps
####

fasting_glucose_PhenoScanner_GWAS <- read.delim("./data/fasting_glucose_PhenoScanner_GWAS.tsv")

fasting_glucose <- fasting_glucose_PhenoScanner_GWAS[,c("rsid","hg19_coordinates","hg38_coordinates","a1", "a2","trait","study", "year", "ancestry", "beta", "se", "p", "direction","n")]
fasting_glucose <- fasting_glucose[fasting_glucose$p<5*10^{-8},]
fasting_glucose <- fasting_glucose[fasting_glucose$study=="MAGIC",]
fasting_glucose <- fasting_glucose[fasting_glucose$year==2012,]
fasting_glucose <- fasting_glucose[fasting_glucose$trait=="Fasting glucose",]

write.csv(fasting_glucose, file="./csv files/fasting_glucose_assoc.csv")


#####
# add to snp list for cardio
####

fasting_glucose<-separate(fasting_glucose,"hg19_coordinates", sep=":", into=c("chr","pos"))
fasting_glucose$chr<-gsub("chr", "", fasting_glucose$chr)
qctoolify(fasting_glucose,chr=chr,pos=pos,"./data/fasting_glucose_snps.txt")

############################################################################################################################

# fasting insulin

fasting_insulin_PhenoScanner_GWAS <- read.delim("./data/fasting_insulin_PhenoScanner_GWAS.tsv")

fasting_insulin <- fasting_insulin_PhenoScanner_GWAS[,c("rsid","hg19_coordinates","hg38_coordinates","a1", "a2","trait","study", "year", "ancestry", "beta","se", "p", "direction","n")]
fasting_insulin <- fasting_insulin[fasting_insulin$p<5*10^{-8},]
fasting_insulin <- fasting_insulin[fasting_insulin$study=="MAGIC",]
fasting_insulin <- fasting_insulin[fasting_insulin$year==2012,]
fasting_insulin <- fasting_insulin[fasting_insulin$trait=="log Fasting insulin",]

write.csv(fasting_insulin, file="./csv files/fasting_insulin_assoc.csv")

# create snp list
fasting_insulin<-separate(fasting_insulin,"hg19_coordinates", sep=":", into=c("chr","pos"))
fasting_insulin$chr<-gsub("chr", "", fasting_insulin$chr)
qctoolify(fasting_insulin,chr=chr,pos=pos,"./data/fasting_insulin_snps.txt")

#############################################################################################################################

# Leptin


leptin_PhenoScanner_GWAS <- read.delim("./data/leptin_PhenoScanner_GWAS.tsv")

leptin <- leptin_PhenoScanner_GWAS[,c("rsid","hg19_coordinates","hg38_coordinates","a1", "a2","trait","study", "year", "ancestry", "beta","se", "p", "direction","n")]
leptin <- leptin[leptin$p<5*10^{-8},]
#### NO GLOBALLY SIGNIFICANT VALUES; find data elsewhere

rm(list=c("leptin_PhenoScanner_GWAS", "leptin"))

leptin_assoc <- read.csv("./data/leptin_assoc.csv")
qctoolify(leptin_assoc,chr=chr,pos=pos,"./data/leptin_snps.txt")

################################################################################################
# HOMA(IR)  -1 SNP only

HOMA <- read.csv("./data/HOMA_IR_assoc.csv")
qctoolify(HOMA,chr=chr,pos=pos,"./data/HOMA_snps.txt")

#####################################################################################################################################
# Imflammatory markers





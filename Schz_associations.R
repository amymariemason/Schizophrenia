######
# Author: Jess Rees
# Date: May 2018
# Aim: obtain the genetic associations with depression outcomes from UK Biobank  
######

rm(list=ls())
setwd("/Users/jessica/Dropbox/201805 - cvd and depression")

library(dplyr)

#################################################
#Read in data 
#################################################

#Read in outcome data for depression, genetic variants and PC data.
depression<-read.csv("most.csv",header = TRUE)
variants<-read.csv("depress_all.csv",header = TRUE)
PCs<-read.table("ukbb_eur_unrel_depress_without.sample",header=TRUE) #what is this file?
PCs<-PCs[-1,]

#Genetic associations from SNPTEST for any 
snptest_any<-read.table("any.out",header=TRUE)


#################################################
#Genetic associations with depression
#################################################

#Get rid of indiviudals with no PC data and convert to numeric variables
print(paste0("There are ",length(which(is.na(PCs$PC1)))," with no PC data"))
PCs<-PCs[!is.na(PCs$PC1),]
PC<-c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCs[,PC]<-apply(PCs[,PC],2,function(x) as.numeric(as.character(x)))

#Merge PC data with genetic data
variants<-merge(variants, PCs, by.x="X", by.y="ID_1")

#Drop non-european and then merge genetic data and depression data
depression_euro<-depression[which(depression$euro==1),]
gen_dep<-merge(depression_euro,variants,by.x="UKB_sample_ID",by.y="X")

#Genetic variants we are interested in and the outcomes 
snps<-colnames(gen_dep)[grep("rs",colnames(gen_dep))]
outcome<-c("mod.x","sev.x", "mdo.x", "any.x","father","mother")

#Obtain the genetic associations for all of the genetic variants for the different depression outcomes. 
#These have been adjusted for the first 10 PCs.  
row_num=seq(0,length(snps)*length(outcome),by=length(snps))
coef=as.data.frame(matrix(data=NA,nrow=length(snps)*length(outcome),ncol=5))
for (j in 1:length(outcome)){
  for (i in 1:length(snps)){
    a<-summary(glm(gen_dep[,outcome[j]]~gen_dep[,snps[i]]+gen_dep$PC1+gen_dep$PC2+gen_dep$PC3+gen_dep$PC4
                   +gen_dep$PC5+gen_dep$PC6+gen_dep$PC7+gen_dep$PC8+gen_dep$PC9+gen_dep$PC10, family = binomial))
    coef[i+row_num[j],1]<-as.character(snps[i])
    coef[i+row_num[j],2]<-outcome[j]
    coef[i+row_num[j],3]<-a$coef[2]
    coef[i+row_num[j],4]<-a$coef[2,2]
    coef[i+row_num[j],5]<-a$coef[2,4]
  }
  print(outcome[j])
}
colnames(coef)<-c("snp","outcome","beta","se","p.value")
coef$outcome<-as.factor(coef$outcome)
levels(coef$outcome)<-c("any","father","mdo","mod","mother","sev")
write.csv(coef,"outcome_associations.csv",row.names = FALSE)


######################################################################################################
#Check the above corresponds with the estimates obtained from SNPTEST 
######################################################################################################

#Consider for the outcome 'any' 
any<-coef[which(coef$outcome=="any.x"),]
any<-merge(any, snptest_any[,c("rsid","frequentist_add_beta_1")],by.x="snp",by.y="rsid")
any$diff<-any$beta-any$frequentist_add_beta_1



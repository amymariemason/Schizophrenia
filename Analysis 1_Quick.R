# Auther = Amy Mason
# Date = Nov 2018
# Purpose = Quick LM tests to examine data
#
###################################################
library(MendelianRandomization)
library(epiDisplay)
setwd("//me-filer1/home$/am2609/My Documents/Programs/CVD_Schizophrenia")

########## analysis of family history of CHD on Schitzophrenia  ##################

## Read family history data
most <- read.csv("~/CVD_depression/snps/most.csv")

## euro: european ancestry  (inclusion criterion)
# Restrict to europeans.
most.euro <- most[most[, 12] == 1, ]
dim(most.euro)

### read schitzophrenia data

schitz <- read.csv("~/CVD_Schizophrenia/inputs/schitz.sample", sep="")

##combine data

most.euro2<-merge(most.euro, schitz[,c("ID_1","sch_1","sch_2"), by.x="UKB_sample_ID", by.y="ID_1"])

## Primary observational analysis of family history of CHD on any level of depression.
# Risk factor.
family <- as.numeric(most.euro2[, "father"] == 1 | most.euro2[, "mother"] == 1)
family2 <- as.numeric(most.euro2[, "father"] + most.euro2[, "mother"])

# Outcome.
sch_basic <- as.numeric(as.character(most.euro2$sch_1))
sch_all<-as.numeric(as.character(most.euro2$sch_2))

# Table

table(family, sch_basic)
summary(table(family, sch_basic))

## Observational association.
fit_basic <- glm(sch_basic~ family , family = "binomial")
summary(fit_basic)
logistic.display(fit_basic, decimal =4)

fit_all <- glm(sch_all~ family , family = "binomial")
summary(fit_all)

#
fit_basic2 <- glm(sch_basic~ family2 , family = "binomial")
summary(fit_basic2)
## Family history no significant effect OR(0.96 (0.83-1.10))

################################ genetic predisposition to CHD as predictor
# Risk factor
Genetic<-as.numeric(most.euro2$COMPOSITE)
# summary in yes/no schitzophrenia
summary(most.euro2[most.euro2$sch_1==1,]$COMPOSITE)
summary(most.euro2[most.euro2$sch_1==0,]$COMPOSITE)



fit_genetic<-glm(sch_basic~Genetic, family="binomial")
summary(fit_genetic)
logistic.display(fit_genetic, decimal =4)

# genetic prediction is not significant: coef -0.002 (-0.008,0.004)


############################################################




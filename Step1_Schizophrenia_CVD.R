#
# Written by Amy Mason
# Part of investigation into effect of history of CVD 
# and other co-morbidity factors on schitzophrenia

#This file includes:
#Mendelian Randomization of BMI on schitzophrenia
#
#
#


##########   Mendelian Randomization of BMI on Depression   ##########

library(MendelianRandomization)

## Set working directory.

## Read some data.

most <- read.csv("most.csv")
## mdo: mmoderate depression only
## any: either moderate or severe depression
## euro: european ancestry  (inclusion criterion)

## Restrict to europeans.
most.euro <- most[most[, 12] == 1, ]
dim(most.euro)


## Primary observational analysis of family history of CHD on any level of depression.

## Risk factor.
family <- as.numeric(most.euro[, 9] == 1 | most.euro[, 10] == 1)

## Outcome.
depression <- most.euro$any

## Observational association.
fit1 <- glm(depression ~ family - 1, family = "binomial")
summary(fit1)
## Family history has a significant observational effect.

########################################

## Now move on to analyze BMI effect on depression.

## SNPs provided to phenoscanner (taken from file "locke_snps_to_qctool").
## We copy-paste the names here for simplicity.

snp.names <- c("rs1558902", "rs6567160", "rs13021737", "rs10938397", "rs543874", "rs2207139", "rs11030104",  "rs3101336", 
"rs7138803", "rs10182181", "rs3888190", "rs1516725", "rs12446632", "rs2287019", "rs16951275", "rs3817334", 
"rs2112347", "rs12566985", "rs3810291", "rs7141420", "rs13078960", "rs10968576", "rs17024393", "rs657452", 
"rs12429545", "rs12286929", "rs13107325", "rs11165643", "rs7903146", "rs10132280", "rs17405819", "rs6091540", 
"rs1016287", "rs4256980", "rs17094222", "rs12401738", "rs7599312", "rs2365389", "rs205262", "rs2820292", 
"rs12885454", "rs9641123", "rs12016871", "rs16851483", "rs1167827", "rs758747", "rs1928295", "rs9925964", 
"rs11126666", "rs2650492", "rs6804842", "rs12940622", "rs7164727", "rs11847697", "rs4740619", "rs492400", 
"rs13191362", "rs3736485", "rs17001654", "rs11191560", "rs2080454", "rs7715256", "rs2176040", "rs1528435", 
"rs2075650", "rs1000940", "rs2033529", "rs11583200", "rs7239883", "rs2836754", "rs9400239", "rs10733682", 
"rs11688816", "rs11057405", "rs9914578", "rs977747", "rs2121279", "rs29941", "rs11727676", "rs3849570", 
"rs9374842", "rs6477694", "rs4787491", "rs1441264", "rs7899106", "rs2176598", "rs2245368", "rs17203016", 
"rs17724992", "rs7243357", "rs16907751", "rs1808579", "rs13201877", "rs2033732", "rs9540493", "rs1460676", "rs6465468")

## Read data from phenoscanner (they were provided in two sets of 50 SNPs and saved in the files below).
bmisnps <- rbind(read.csv("BMISNPs.Phenoscanner.csv"), read.csv("BMISNPs2.Phenoscanner.csv"))

## Restrict the data matrix to BMI associations and europeans.
bmisnps.to.use <- bmisnps[bmisnps$Trait == "BMI" & bmisnps$Ancestry == "European", ]
dim(bmisnps.to.use)   ## Twice the number of SNPs.

bmisnps.to.use$N.Studies   ## Alternates between 46 and 114. Better to restrict to the biggest studies.

bmisnps.to.use <- bmisnps.to.use[bmisnps.to.use$N.Studies == 114, ]
dim(bmisnps.to.use)   ## half.


## Cross-check. that everything ran smoothly and that SNPs are in the same order as in snp.names.
bmisnps.to.use$rsID %in% snp.names
bmisnps.to.use$rsID == snp.names
## One SNP (rs9581854) did not appear in the original list (it is replaced by rs12016871).
## Possibly the two SNPs are in LD and phenoscanner used the first as a proxy for the second?

## Update SNP names anyway:
snp.names <- as.character(bmisnps.to.use$rsID)

## Take snp-risk factor effects, standard errors, mafs and effect alleles from the phenoscanner runs.
bx <-bmisnps.to.use$Beta
sx <- bmisnps.to.use$SE
maf.x <- bmisnps.to.use$MAF
effect.allele.x <- bmisnps.to.use$Effect.Allele

########################################

## Now on to depression.

## Start with any type of depression. Get data from any.out.
any.out <- read.table("any.out", header = TRUE, sep = " ")

## Restrict to the SNPs studying BMI.
any.out <- any.out[any.out$rsid %in% snp.names, ]
dim(any.out)   ## 97 rows, as many as the SNPs used for BMI.

outcome.snp.rsid <- as.character(any.out$rsid)
outcome.snp.rsid %in% snp.names   ## all TRUE. But the SNPs are not in the same order.

## To fix this, do
match.snps <- match(snp.names, outcome.snp.rsid)
snp.names == outcome.snp.rsid[match.snps]   ## fine.

## This re-orders any.out to match the order of SNPs in "snp.names".
any.out <- any.out[match.snps, ]
## To verify:
any.out$rsid[1:10]
bmisnps.to.use$rsID[1:10]

## Get the SNP-outcome associations from "any.out".
by <- any.out$frequentist_add_beta_1
sy <- any.out$frequentist_add_se_1
py <- any.out$frequentist_add_pvalue
maf.y <- any.out$all_maf
effect.allele.y <- any.out$alleleB

## Summarize results.
results <- cbind(snp.names, bx, sx, maf.x, as.character(effect.allele.x), by, sy, py, maf.y, as.character(effect.allele.y))
colnames(results) <- c("rsid", "bx", "sx", "maf.x", "eff.x", "by", "sy", "py", "maf.y", "eff.y")



## When effect alleles are different, need to flip the sign of the corresponding effect estimate.

## Cross-check that the effect allele in the BMI dataset is always one of the two alleles in any.out.
( as.character(effect.allele.x) == as.character(effect.allele.y) ) | ( as.character(effect.allele.x) == as.character(any.out$alleleA) )

## Compare the effect alleles to decide when to flip signs, then do the flipping.
when.to.flip <- as.character(effect.allele.x) != as.character(effect.allele.y)
by[when.to.flip] <- - by[when.to.flip]
## Note: we flipped only the vector by, not the results table or any.out.

########################################

## Now do Mendelian randomization.

mr.fit <- mr_input(bx, sx, by, sy)

mr_ivw(mr.fit)
mr_allmethods(mr.fit, method = "main")

## There seems to be no significant effect of BMI on (any level of) depression.

########################################

## We now repeat the same procedure for moderate only (mdo) and severe only (sev~) depression.

## Get data from the corresponding files.
mdo.out <- read.table("mdo.out.txt", header = TRUE, sep = " ")
sev.out <- read.table("sev.out.txt", header = TRUE, sep = " ")

## Restrict to the SNPs studying BMI.
mdo.out <- mdo.out[mdo.out$rsid %in% snp.names, ]
sev.out <- sev.out[sev.out$rsid %in% snp.names, ]

## Cross-check that these have the same SNPs and in the same order as the original any.out.
outcome.snp.rsid == mdo.out$rsid
outcome.snp.rsid == sev.out$rsid

## They do, so apply the same reordering to them as we did for any.out.
mdo.out <- mdo.out[match.snps, ]
sev.out <- sev.out[match.snps, ]

## Get the SNP-outcome associations for moderate only depression.
by.mdo <- mdo.out$frequentist_add_beta_1
sy.mdo <- mdo.out$frequentist_add_se_1
py.mdo <- mdo.out$frequentist_add_pvalue
## Minor allele frequencies and effect alleles are the same.

results.mdo <- cbind(snp.names, bx, sx, maf.x, as.character(effect.allele.x), by.mdo, sy.mdo, py.mdo, maf.y, as.character(effect.allele.y))
colnames(results.mdo) <- c("rsid", "bx", "sx", "maf.x", "eff.x", "by", "sy", "py", "maf.y", "eff.y")


## Likewise for severe only depression.
by.sev <- sev.out$frequentist_add_beta_1
sy.sev <- sev.out$frequentist_add_se_1
py.sev <- sev.out$frequentist_add_pvalue

results.sev <- cbind(snp.names, bx, sx, maf.x, as.character(effect.allele.x), by.sev, sy.sev, py.sev, maf.y, as.character(effect.allele.y))
colnames(results.sev) <- c("rsid", "bx", "sx", "maf.x", "eff.x", "by", "sy", "py", "maf.y", "eff.y")

## Flip the signs of the SNP-outcome effects when effect alleles differ.
by.mdo[when.to.flip] <- - by.mdo[when.to.flip]
by.sev[when.to.flip] <- - by.sev[when.to.flip]

########################################

## Now do Mendelian randomization for moderate depression only.

mr.fit.mdo <- mr_input(bx, sx, by.mdo, sy.mdo)
mr_ivw(mr.fit.mdo)
mr_allmethods(mr.fit.mdo, method = "main")
## There seems to be no significant effect of BMI on moderate depression.

## Finally, do Mendelian randomization for severe depression only.

mr.fit.sev <- mr_input(bx, sx, by.sev, sy.sev)
mr_ivw(mr.fit.sev)
mr_allmethods(mr.fit.sev, method = "main")
## There seems to be no significant effect of BMI on severe depression.

########################################



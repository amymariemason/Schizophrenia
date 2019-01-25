# Auther = Amy Mason
# Date = Nov 2018
# Purpose = create snp list for using on cardio
#
###################################################

library(readxl)
setwd("//me-filer1/home$/am2609/My Documents")

##################################################
# load in each set of snps for: imflammatory, lipids, bmi, bp and insulin types
# keep all of the chr/pos in running set called "all_SNPs"


# Inflammatory snps
inflamatory<- read_excel("~/CVD_Schizophrenia/data/inflamatory_markers_v2.xlsx")
inflamatory <- cbind(inflamatory, do.call("rbind", strsplit(inflamatory$`chr:locus`, ":")))
names(inflamatory)[11:12]<-c("chr","pos")
names(inflamatory)[4]<-"rsid"
inflamatory$chr<-gsub("chr","",inflamatory$chr)


all_SNPs<-unique(inflamatory[,c("rsid", "chr", "pos")])

# lipids

lipids<- read.csv("~/CVD_Schizophrenia/data/lipid_mrcat_all_MR_Catalogue_v2.csv")
lipids<-lipids[lipids$Source=="GLGC",]
lipids<-lipids[lipids$Year=="2010",]
names(lipids)[2]<-"rsid"
names(lipids)[3]<-"pos_name"

lipids <- cbind(lipids, do.call("rbind", strsplit(as.character(lipids$pos_name), ":")))
names(lipids)[25]<-"chr"
names(lipids)[26]<-"pos"
lipids$chr<-gsub("chr", "", lipids$chr)

all_SNPs<-rbind(all_SNPs, lipids[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)


###bmi snps
# this list came from Stephen so matches Hackathon data

bmi_snps = scan(what="character")
rs1558902
rs6567160
rs13021737
rs10938397
rs543874
rs2207139
rs11030104
rs3101336
rs7138803
rs10182181
rs3888190
rs1516725
rs12446632
rs2287019
rs16951275
rs3817334
rs2112347
rs12566985
rs3810291
rs7141420
rs13078960
rs10968576
rs17024393
rs657452
rs12429545
rs12286929
rs13107325
rs11165643
rs7903146
rs10132280
rs17405819
rs6091540
rs1016287
rs4256980
rs17094222
rs12401738
rs7599312
rs2365389
rs205262
rs2820292
rs12885454
rs9641123
rs12016871
rs16851483
rs1167827
rs758747
rs1928295
rs9925964
rs11126666
rs2650492
rs6804842
rs12940622
rs7164727
rs11847697
rs4740619
rs492400
rs13191362
rs3736485
rs17001654
rs11191560
rs2080454
rs7715256
rs2176040
rs1528435
rs2075650
rs1000940
rs2033529
rs11583200
rs7239883
rs2836754
rs9400239
rs10733682
rs11688816
rs11057405
rs9914578
rs977747
rs2121279
rs29941
rs11727676
rs3849570
rs9374842
rs6477694
rs4787491
rs1441264
rs7899106
rs2176598
rs2245368
rs17203016
rs17724992
rs7243357
rs16907751
rs1808579
rs13201877
rs2033732
rs9540493
rs1460676
rs6465468



bmi= read.csv("~/CVD_Schizophrenia/data/lipid_mrcat_all_MR_Catalogue_v2.csv")
bmi<-bmi[bmi$PMID=="25673413",]
bmi<-bmi[bmi$Trait=="BMI",]
# only 2 of the above listed SNPS are in here, so use the bmi set you have data
# rather than the list above
#bmi<-bmi[bmi$rsID%in%bmi_snps,]

#phenoscanner data for snp list for bmi
bmi2<-read.delim("~/CVD_Schizophrenia/data/bmi_rs_PhenoScanner.tsv")
names(bmi2)[5]<-"Effect.allele"
names(bmi2)[6]<-"Association.alleles"
bmi2<-bmi2[bmi2$pmid=="25673413",]
bmi2<-bmi2[bmi2$ancestry=="European",]
bmi2<-bmi2[bmi2$trait=="Body mass index",]
## this is now the same set as Steve's so use this. 


bmi2 <- cbind(bmi2, do.call("rbind", strsplit(as.character(bmi2$hg19_coordinates), ":")))
names(bmi2)[23]<-"chr"
names(bmi2)[24]<-"pos"
bmi2$chr<-gsub("chr", "", bmi2$chr)

all_SNPs<-rbind(all_SNPs, bmi2[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)




#blood pressure

bp <- read.csv("~/CVD_Schizophrenia/data/bp_associations.csv")
names(bp)[1]="chr"
all_SNPs$pos<-as.character((all_SNPs$pos))
all_SNPs<-rbind(all_SNPs, bp[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)



# fasting glucose
fasting_glucose <- read.csv("~/CVD_Schizophrenia/data/fasting_glucose_assoc.csv")
names(fasting_glucose)[5]="Effect.allele"
names(fasting_glucose)[6]="Association.allele"

fasting_glucose <- cbind(fasting_glucose, do.call("rbind", strsplit(as.character(fasting_glucose$hg19_coordinates), ":")))
names(fasting_glucose)[16]<-"chr"
names(fasting_glucose)[17]<-"pos"
fasting_glucose$chr<-gsub("chr", "", fasting_glucose$chr)

all_SNPs<-rbind(all_SNPs, fasting_glucose[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)

# fasting insulin

fasting_insulin <- read.csv("~/CVD_Schizophrenia/data/fasting_insulin_assoc.csv")
names(fasting_insulin)[5]="Effect.allele"
names(fasting_insulin)[6]="Association.allele"
miss
fasting_insulin <- cbind(fasting_insulin, do.call("rbind", strsplit(as.character(fasting_insulin$hg19_coordinates), ":")))
names(fasting_insulin)[16]<-"chr"
names(fasting_insulin)[17]<-"pos"
fasting_insulin$chr<-gsub("chr", "", fasting_insulin$chr)

all_SNPs<-rbind(all_SNPs, fasting_insulin[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)

# HOMA- IR

HOMA_IR <- read.csv("~/CVD_Schizophrenia/data/HOMA_IR_assoc.csv")
names(HOMA_IR)[2]="Effect.allele"
names(HOMA_IR)[3]="Association.allele"
HOMA_IR$Association.allele<-"T"

all_SNPs<-rbind(all_SNPs, HOMA_IR[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)

#leptin

leptin<- read.csv("~/CVD_Schizophrenia/data/leptin_assoc.csv")

all_SNPs<-rbind(all_SNPs, leptin[,c("rsid", "chr", "pos")])
all_SNPs<-unique(all_SNPs)

#########################################################
# output : create text document of snps suitable for cardio



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


qctoolify(all_SNPs,filename="~/CVD_Schizophrenia/inputs/all_snps.txt")



##################### repeat of extra snps for fibrinogen

library(readxl)
setwd("//me-filer1/home$/am2609/My Documents")


# Inflammatory snps
Fibrinogen_PhenoScanner_GWAS <- read.delim("~/CVD_Schizophrenia/data/Fibrinogen_PhenoScanner_GWAS.tsv")
Fibri <-Fibrinogen_PhenoScanner_GWAS[!is.na(Fibrinogen_PhenoScanner_GWAS$beta),] 
Fibri <-Fibri[Fibri$ancestry=="European",] 
Fibri <-Fibri[Fibri$p<5*10^(-8),] 
Fibri <- cbind(Fibri, do.call("rbind", strsplit(as.character(Fibri$hg19_coordinates), ":")))
names(Fibri)[32:33]<-c("chr","pos")
Fibri$chr<-gsub("chr","",Fibri$chr)


all_SNPs2<-unique(Fibri[,c("rsid", "chr", "pos")])
qctoolify(all_SNPs2,filename="~/CVD_Schizophrenia/inputs/extra__fibro_snps.txt")

###################### all in one for PGC

all_SNPs<-rbind(all_SNPs, Fibri[,c("rsid", "chr", "pos")])

# load PGC snps (downloaded from https://www.med.unc.edu/pgc/results-and-downloads/ , SCZ2, all SNPs on 18th Jan 2019)
rall <- read.delim("~/CVD_Schizophrenia/PGC/rall.txt")
PGC_snps<-rall[rall$snpid %in% all_SNPs$rsid,]

write.csv(PGC_snps, file="~/CVD_Schizophrenia/inputs/PGC_cutdown.csv")

# 6 snps missing:  rs1998013, rs2247056, rs11246602, rs139385870, rs36083386, rs11442819  
# find proxies?
proxy_PGC <- read.csv("~/CVD_Schizophrenia/proxy_PGC.csv")
PGC_snps2<- rall[rall$snpid %in% proxy_PGC$Proxy.snp,]

write.csv(PGC_snps2, file="~/CVD_Schizophrenia/inputs/PGC_cutdown_proxy.csv")

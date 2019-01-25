# Auther = Amy Mason
# Date = Nov 2018
# Purpose = create Schitzophrenia outcomes sample set for using on cardio
#
###################################################


# load adiposity sample set
setwd("//me-filer1/home$/am2609/My Documents/Programs/GWAS_inprogress/Outputs")
sample_all<- read.csv("~/Programs/GWAS_inprogress/Outputs/AdiposityID_sampleset", sep="")



# load schitzophrenia outcomes

outcomes <- read.csv("~/CVD_Schizophrenia/inputs/schitzophrenia_outcomes.csv")
wanted_outcomes =c("sch_1", "sch_2")

#########
# create function that returns column of binary outcomes for whole 
#empty results files
output<-rep(NA, nrow(sample_all))
output<-as.data.frame(output)

extract_outcome<-function(sample_list, sample_id="ID_1", outcomefile, outcome, outcome_id="eid"){
  stopifnot(!is.null(outcomefile[,outcome]))
  stopifnot(!is.null(outcomefile[,outcome_id]))
  stopifnot(!is.null(sample_list[,sample_id]))
  events<-which(outcomefile[,outcome]==1)
  cases = unique(outcomefile[events,outcome_id])
  output_temp = ifelse(as.numeric(sample_list[2:nrow(sample_list),sample_id])%in%cases, 1, 0)
  output_temp<-as.data.frame(output_temp)
  output_temp<-rbind("B", output_temp)
  names(output_temp)<-outcome
  return(output_temp)
}

# Apply this accross all 

allframes = lapply(wanted_outcomes,
                   function(x)extract_outcome(sample_list=sample_all, 
                                              sample_id="ID_1", 
                                              outcomefile=outcomes, 
                                              outcome=x, 
                                              outcome_id="n_eid"))
answer = do.call(cbind,allframes)

# check same row number and bind
assertthat::are_equal(nrow(sample_all),nrow(answer))
sample_out<- cbind(sample_all, answer)



#return to using UKBB ID's
sample_out[,1]<-sample_out[,2]

#check row order is same as original sample file
sample_out$id2<-1:nrow(sample_out)
assertthat::are_equal(sample_out$id, sample_out$id2)

#remove unneeded columns
sample_out<-sample_out[, !names(sample_out) %in% c("id", "exclude", "error_check","id2")]

# fix error with snptest not reading missing correctly
sample_out[1, "missing"]<-0

# save sample file
write.table(sample_out, "~/CVD_Schizophrenia/inputs/schitz.sample", row.names=FALSE, quote=FALSE)



library(rjags)
library(postjags)
library(mcmcplots)
setwd("C:\\Users\\mkf58\\Google Drive\\Oglelab_NAU_JOB\\NAU_JOB\\MultiComp\\Review_Example\\MulticompExample\\Guo_Jornada\\Monsoon")

updateits <- 1000

hbits <- 10000
hbthin <- 1

cpits <- 10000
cpthin <- 1

#nhits <- 10000
#nhthin <- 1

nhits <- 40000
nhthin <- 20


#hbits <- 100
#hbthin <- 1

#cpits <- 100
#cpthin <- 1

#nhits <- 100
#nhthin <- 1

#load("inits_baseline.r")
load("inits_baseline_nh.Rdata")

load.module('dic')
#setwd("Y:/common/Projects/mult_comparison/RealExamples/Guo_Jornada")

wp_all<-read.csv("Data_wp_all.csv", header=T)
str(wp_all)
wp_all$dt<-as.POSIXct(paste(wp_all$dt), format="%m/%d/%Y %H:%M")
baseline<-subset(wp_all, pd==1)
#extract only the odd-numbered shrubs with sufficient sample size
baseline<-subset(baseline, shrubID==1|shrubID==3|shrubID==5|shrubID==7|shrubID==9|shrubID==11|shrubID==13|shrubID==15)
baseline$shrubID2<-c()
for(i in 1:nrow(baseline)){
  baseline$shrubID2[i]<-if(baseline$shrubID[i]==1){1} else
    if(baseline$shrubID[i]==3){2} else
      if(baseline$shrubID[i]==5){3} else
        if(baseline$shrubID[i]==7){4} else
          if(baseline$shrubID[i]==9){5} else
            if(baseline$shrubID[i]==11){6}else
              if(baseline$shrubID[i]==13){7} else {8}
}
attach(baseline)
str(baseline)
#lowest value
min(baseline$wp)

#set means for covariate centering of regression parameters
meanmaxD<-5
meanswc30<-10
meanswc60<-13

#number of observations
N<-nrow(baseline)
#for delta trick for all antecedent parameters
alphaA <- c(1,1,1,1,1,1,1)
alphaB <- c(1,1,1,1,1,1,1)
alphaC <- c(1,1,1,1,1,1,1)

#create list of data objects for model;
#Note: extracts diurnal data from the wp_all data (which also includes predawn or baseline data)
datalist<-list(wp=wp, shrubID=shrubID2,
               maxD=as.matrix(baseline[,match("maxD1",colnames(baseline)):match("maxD7", colnames(baseline))]), 
               swc30=as.matrix(baseline[,match("swc30_1",colnames(baseline)):match("swc30_7", colnames(baseline))]), 
               swc60=as.matrix(baseline[,match("swc60_1",colnames(baseline)):match("swc60_7", colnames(baseline))]),
               meanmaxD=meanmaxD, meanswc30=meanswc30, 
               meanswc60=meanswc60,
               N=N, alphaA=alphaA, alphaB=alphaB, alphaC=alphaC,
               Ab=rep(10,7))

#create initials
# inits<-list(list(), list(), list())
# for(i in 1:3){
#   inits[[i]][[1]]<-saved.state.nh[[2]][[i]][[1]][,c(1,3,5,7,9,11,13,15)]
#   inits[[i]][[2]]<-saved.state.nh[[2]][[i]][[2]]
#   inits[[i]][[3]]<-saved.state.nh[[2]][[i]][[3]]
#   inits[[i]][[4]]<-saved.state.nh[[2]][[i]][[4]]
#   inits[[i]][[5]]<-saved.state.nh[[2]][[i]][[5]]
#   names(inits[[i]])<-c("b", "deltaA", "deltaB", "deltaC", "tau")
# }

# Load the initial values
load("inits_baseline_nh.Rdata")
### For all models, need to monitor log-pointwise density (lpd), to compute PPI. Currently not
### monitored because lpd is a large vector (number of observations)

### Run NH model
jm_nh <- jags.model("R model baseline NH.R", data=datalist, inits=saved.state.nh[[2]], n.chains=3, n.adapt = updateits)

runtime.nh <- system.time(
  {
    #run model
    coda_nh <- coda.samples(jm_nh, variable.names = c("deviance", "tau", "sig",
                                                      "wA", "wB", "wC", "log_lik",
                                                      "b", "deltaA", "deltaB", "deltaC"),
                            n.iter=nhits, thin=nhthin)
  }
)
#save model
save(jm_nh, file = "Output/jm_nh_baseline.Rdata")
save(coda_nh, file = "Output/coda_nh_baseline.Rdata")
print(runtime.nh)
mcmcplot(coda_nh, parms=c("deviance", "tau", "sig", "tau.eps.b",
                         "wA", "wB", "wC",
                         "b", "mu.b"))
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(coda_nh[[1]]))
sel <- c(sel, grep(pattern="log_lik", x=colnames(coda_nh[[1]])))
for(i in 1:3){
    gelcoda[[i]] <- as.mcmc(coda_nh[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

#extract final iteration
newinits<-initfind(coda_nh)#default is final iteration
newinits[[1]]
#remove non-root node variables
#retain a, deltaA-E, tau
saved.state.nh <- removevars(initsin = newinits, variables=c(5:6, 8:10)) # remove non-variable nodes
#check both items in list
#saved.state.nh[[1]]
saved.state.nh[[2]]#this goes into jags.model() to reinitialize
save(saved.state.nh, file="inits_baseline_nh.Rdata") #for local
#load(file="inits_baseline_nh.Rdata")

#check convergence
gel.nh<-gelman.diag(coda_nh, multivariate = F)
gel.nh$psrf[match("deviance", row.names(gel.nh$psrf)),]
#by param
b_params<-matrix(NA, nrow=8, ncol=7)
for(i in 1:8){
  for(j in 1:7){
    b_params[i,j]<-paste0("b[", j, ",", i, "]")
  }
}
#gel.nh$psrf[match(b_params[,1], row.names(gel.nh$psrf)),]
#gel.nh$psrf[match(b_params[,2], row.names(gel.nh$psrf)),]
#gel.nh$psrf[match(b_params[,3], row.names(gel.nh$psrf)),]
#gel.nh$psrf[match(b_params[,4], row.names(gel.nh$psrf)),]
#gel.nh$psrf[match(b_params[,5], row.names(gel.nh$psrf)),]
#gel.nh$psrf[match(b_params[,6], row.names(gel.nh$psrf)),]
#el.nh$psrf[match(b_params[,7], row.names(gel.nh$psrf)),]

#gel.nh$psrf[match("wA[1]", row.names(gel.nh$psrf)):match("wA[7]", row.names(gel.nh$psrf)),]
#gel.nh$psrf[match("wB[1]", row.names(gel.nh$psrf)):match("wB[7]", row.names(gel.nh$psrf)),]
#gel.nh$psrf[match("wC[1]", row.names(gel.nh$psrf)):match("wC[7]", row.names(gel.nh$psrf)),]
#gel.nh$psrf[match("tau", row.names(gel.nh$psrf)),]

save.image(file='Monsoonworkspace_nh.RData')
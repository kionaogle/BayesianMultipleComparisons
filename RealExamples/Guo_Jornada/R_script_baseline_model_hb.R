library(rjags)
library(postjags)
library(mcmcplots)
setwd("C:\\Users\\mkf58\\Google Drive\\Oglelab_NAU_JOB\\NAU_JOB\\MultiComp\\Review_Example\\MulticompExample\\Guo_Jornada\\Monsoon")

updateits <- 1000
hbits <- 2000000
hbthin <- 1000

#cpits <- 10000
#cpthin <- 1

#nhits <- 10000
#nhthin <- 1

#hbits <- 100
#hbthin <- 1

#cpits <- 100
#cpthin <- 1

#nhits <- 100
#nhthin <- 1

#load("inits_baseline.r")
load("inits_baseline_hb.Rdata")

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


### For all models, need to monitor log-pointwise density (lpd), to compute PPI. Currently not
### monitored because lpd is a large vector (number of observations)

### Run NH model
runtime.hb <- system.time(
  {
    jm_hb <- jags.model("R model baseline HB.R", data=datalist, inits=saved.state.hb[[2]], n.chains=3, n.adapt = updateits)
    #run model
    coda_hb <- coda.samples(jm_hb, variable.names = c("deviance", "tau", "sig", "tau.eps.b",
                                                      "wA", "wB", "wC", "lpd", "b.control", "b.watered", "b.diff",
                                                      "b", "mu.b", "deltaA", "deltaB", "deltaC", "p.b"),
                            n.iter=hbits, thin=hbthin)
  }
)
#save model
save(jm_hb, file = "Output/jm_hb_baseline.Rdata")
save(coda_hb, file = "Output/coda_hb_baseline.Rdata")
print(runtime.hb)
mcmcplot(coda_hb, parms=c("deviance", "tau", "sig", "tau.eps.b",
                         "wA", "wB", "wC",
                         "b", "mu.b"))
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(coda_hb[[1]]))
sel <- c(sel, grep(pattern="log_lik", x=colnames(coda_hb[[1]])))
for(i in 1:3){
    gelcoda[[i]] <- as.mcmc(coda_hb[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

hbsum <- coda.fast(chains=3, coda=coda_hb)
#extract final iteration
newinits<-initfind(coda_hb)#default is final iteration
newinits[[1]]
#remove non-root node variables
#retain a, deltaA-E, tau
saved.state.hb <- removevars(initsin = newinits, variables=c(1,5,7,10:12 )) # remove non-variable nodes
#check both items in list
#saved.state.hb[[1]]
saved.state.hb[[2]]#this goes into jags.model() to reinitialize
save(saved.state.hb, file="inits_baseline_hb2.Rdata") #for local
#load(file="inits_baseline_hb.Rdata")

#check convergence
#gel.hb<-gelman.diag(coda_hb, multivariate = F)
#gel.hb$psrf[match("deviance", row.names(gel.hb$psrf)),]
#gel.hb$psrf[match("b[1,1]", row.names(gel.hb$psrf)):match("b[7,16]", row.names(gel.hb$psrf)),]
#gel.hb$psrf[match("wA[1]", row.names(gel.hb$psrf)):match("wA[7]", row.names(gel.hb$psrf)),]
#gel.hb$psrf[match("wB[1]", row.names(gel.hb$psrf)):match("wB[7]", row.names(gel.hb$psrf)),]
#gel.hb$psrf[match("wC[1]", row.names(gel.hb$psrf)):match("wC[7]", row.names(gel.hb$psrf)),]
#gel.hb$psrf[match("tau", row.names(gel.hb$psrf)),]

save.image(file='Monsoonworkspace_hb.RData')
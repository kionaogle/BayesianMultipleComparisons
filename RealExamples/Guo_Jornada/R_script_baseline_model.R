###############################################################################
# Load the necessary libraries
# Also REMEMBER TO SET THE WORKING DIRECTORY, the folder containing 
# this file.
###############################################################################

if("rjags" %in% installed.packages()[,1] == FALSE){
	install.packages("rjags", repos="https://cran.cnr.berkeley.edu/")
	library(rjags)
} else{
	library(rjags)
}
if("coda" %in% installed.packages()[,1] == FALSE){
	install.packages("coda", repos="https://cran.cnr.berkeley.edu/")
	library(coda)
} else{
	library(coda)
}
if("coda" %in% installed.packages()[,1] == FALSE){
	install.packages("coda", repos="https://cran.cnr.berkeley.edu/")
	library(mcmcplots)
} else{
	library(mcmcplots)
}

###############################################################################
# Variables controlling model runs
###############################################################################

# Set variables controlling the number of update iterations to update the model
# when initializing it, the length of the burnin for the MCMC chains, and how
# many iterations to run the model for.
# n.chains - The number of chains in the MCMC
# n.burnin - The length of burnin to be used
# adaptits - The number of iterations for adaptation (tuning) the MCMC
# hb\cp\nhits - the number of iterations the model will be ran
# hb\cp\nhthin - the thinning interval to be applied when running the model.
n.chains <- 3
n.burnin <- 0
adaptits <- 1000

hbn.burnin <- 0
hbits <- 2000000
hbthin <- 1000

cpn.burnin <- 0
cpits <- 20000
cpthin <- 10

nhn.burnin <- 0
nhits <- 40000
nhthin <- 20

###############################################################################
# Hierarchical Model
###############################################################################
# Load initial values 
load("inits/inits_baseline_hb.Rdata")
load.module('dic')

# Load data and reshape it for running the model
wp_all<-read.csv("data/Data_wp_all.csv", header=T)
#str(wp_all)
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
#str(baseline)
#lowest value
#min(baseline$wp)

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

### For all models, need to monitor log-pointwise density (lpd), to compute PPI. 
### Run HB model
runtime.hb <- system.time(
  {
    jm_hb <- jags.model("models/R model baseline HB.R", data=datalist, inits=saved.state.hb[[2]], n.chains=n.chains, n.adapt = adaptits)
    #run model
    coda_hb <- coda.samples(jm_hb, variable.names = c("deviance", "tau", "sig", "tau.eps.b",
                                                      "wA", "wB", "wC", "lpd", "b.control", "b.watered", "b.diff",
                                                      "b", "mu.b", "deltaA", "deltaB", "deltaC", "p.b"),
                            n.iter=hbits, thin=hbthin, n.burnin=hbn.burnin)
  }
)
#save model output
save(coda_hb, file = "output/coda_hb_baseline.Rdata")
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

# uncomment the below line to save the workspace
#save.image(file='Jornada_hb_workspace.RData')

###############################################################################
# Complete Pooling Model
###############################################################################

# Load initial values 
load("inits/inits_baseline_cp.Rdata")
load.module('dic')

wp_all<-read.csv("data/Data_wp_all.csv", header=T)
#str(wp_all)
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
#str(baseline)
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

### For all models, need to monitor log-pointwise density (lpd), to compute PPI. 
### Run CP model               
runtime.cp <- system.time(
  {
    jm_cp <- jags.model("models/R model baseline CP.R", data=datalist, inits=saved.state.cp[[2]], n.chains=n.chains, n.adapt = adaptits)
    #run model
    coda_cp <- coda.samples(jm_cp, variable.names = c("deviance", "tau", "sig",
                                                      "wA", "wB", "wC", "log_lik",
                                                      "b", "deltaA", "deltaB", "deltaC"),
                            n.iter=cpits, thin=cpthin, cpn.burnin)
  }
)

#save model
save(coda_cp, file = "output/coda_cp_baseline.Rdata")
print(runtime.cp)
mcmcplot(coda_cp, parms=c("deviance", "tau", "sig", "tau.eps.b",
                         "wA", "wB", "wC",
                         "b", "mu.b"))
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(coda_cp[[1]]))
sel <- c(sel, grep(pattern="log_lik", x=colnames(coda_cp[[1]])))
for(i in 1:3){
    gelcoda[[i]] <- as.mcmc(coda_cp[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

# uncomment the below line to save the workspace
#save.image(file='Monsoonworkspace_cp.RData')

###############################################################################
# Non-Hierarchical Model
###############################################################################
#load("inits_baseline.r")
load("inits/inits_baseline_nh.Rdata")
load.module('dic')

wp_all<-read.csv("data/Data_wp_all.csv", header=T)
#str(wp_all)
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
#str(baseline)
#lowest value
#min(baseline$wp)

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
               N=N, alphaA=alphaA, alphaB=alphaB, alphaC=alphaC)

### For all models, need to monitor log-pointwise density (lpd), to compute PPI. Currently not
### monitored because lpd is a large vector (number of observations)

### Run NH model
runtime.nh <- system.time(
  {
    jm_nh <- jags.model("models/R model baseline NH.R", data=datalist, inits=saved.state.nh[[2]], n.chains=3, n.adapt = hbn.burnin)

    #run model
    coda_nh <- coda.samples(jm_nh, variable.names = c("deviance", "tau", "sig",
                                                      "wA", "wB", "wC", "log_lik",
                                                      "b", "deltaA", "deltaB", "deltaC"),
                            n.iter=nhits, thin=nhthin, hbn.burnin)
  }
)
#save model
save(coda_nh, file = "output/coda_nh_baseline.Rdata")
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

#check convergence
gel.nh<-gelman.diag(coda_nh, multivariate = F)
gel.nh$psrf[match("deviance", row.names(gel.nh$psrf)),]

# uncomment the below line to save the workspace
#save.image(file='Monsoonworkspace_nh.RData')


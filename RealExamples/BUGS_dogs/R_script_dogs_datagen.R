#setwd("C:\\Users\\mkf58\\Google Drive\\Oglelab_NAU_JOB\\NAU_JOB\\MultiComp\\Review_Example\\MulticompExample\\BUGS_dogs")

# Script for running the dog shock models.
if("rjags" %in% installed.packages()[,1] == FALSE){
	install.packages("rjags")
	library(rjags)
} else{
	library(rjags)
}
if("mcmcplots" %in% installed.packages()[,1] == FALSE){
	install.packages("mcmcplots")
	library(mcmcplots)
} else{
	library(mcmcplots)
}
if("postjags" %in% installed.packages()[,1] == FALSE){
    install.packages("devtools")
    library(devtools)
    devtools::install_github("fellmk/PostJAGS/postjags")
    library(postjags)
} else{
	library(postjags)
}
load.module('dic')

#load in data
d <- read.csv("dogshock_matrix.csv", header = F)
Y = as.matrix(d)
dim(Y)

Ndog=nrow(d)
Ncomp = Ndog*(Ndog-1)/2

# Indices for dog pairwise comparisons; number of comparisons = (Ndog)*(Ndog-1)/2
pairID = matrix(data =NA, nrow = Ncomp, ncol=2)
ii = 0
for(s1 in 1:(Ndog-1)){
  for(s2 in (s1+1):Ndog){
    ii = ii +1
    pairID[ii,1] = s1
    pairID[ii,2] = s2
  }
}

#outside of JAGS, create the xa and xs matrices of covariates (total past successes & total past failures)
xa <- matrix(NA, nrow=30, ncol=25)
xs <- matrix(NA, nrow=30, ncol=25)
for(i in 1:30){
  #set starting values for the past number of successes (avoiding shock)
  xa[i, 1] <- 0
  #for the past number of failures (getting shocked)
  xs[i, 1] <- 0
    for(t in 2:25){
    #sum number of past successes (avoiding shock) prior to trial t
    xa[i, t] <- sum(Y[i, 1:(t-1)])
    #sum number of past failures (shocked) prior to trial t
    xs[i, t] <- (t-1) - xa[i,t]
  }
}

#create data list; set d as matrix
dat1 <- list(Y=as.matrix(d), Ndog=nrow(d), Ntrial=ncol(d), xa=xa, xs=xs, pairID = pairID, Ncomp = Ncomp)
data <- dat1
########################################
### HB model

# initials
initsHB <-list(list(mu.log.alpha = -1, mu.log.beta = -1, tau.log.beta = 10, tau.log.alpha = 10),
                 list(mu.log.alpha = 1, mu.log.beta = 0, tau.log.beta = 20, tau.log.alpha = 50),
                 list(mu.log.alpha = -.1, mu.log.beta = .8, tau.log.beta = 1, tau.log.alpha = 5))


jmHB <-jags.model("models/R_model_dogs_HB.R", data=dat1, inits=initsHB, n.chains = 3)

#rerun coda.samples (will need to monitor lpd)
runtime.hb <- system.time(
{
update(jmHB, 8000)
codaHB<-coda.samples(jmHB, variable.names = c("deviance", "alpha", "beta", "alpha.diff","beta.diff",
                                              "mu.alpha","mu.beta","tau.log.alpha","tau.log.beta",
                                              "p.alpha", "p.beta", "lpd"), n.iter=5000)
}
)
save(codaHB, file = "output/jm_hb_lpd.Rdata")
#load("Output/jm_hb2lpd.Rdata")
#load("Output/jm_hb.Rdata")

########################################
### CP model

#initials
initsCP=list(list(log.alpha = -1, log.beta = -1), 
               list(log.alpha = -2, log.beta = -2),
               list(log.alpha = -3.5, log.beta = -.2))

# initialize model
jmCP <-jags.model("R_model_dogs_CP.R", data=dat1, inits=initsCP, n.chains = 3)

#rerun coda.samples (will need to monitor lpd)
runtime.cp <- system.time(
{
update(jmCP, 8000)
codaCP<-coda.samples(jmCP, variable.names = c("deviance", "alpha", "beta", "alpha.diff","beta.diff",
                                               "mu.alpha","mu.beta","tau.log.alpha","tau.log.beta", "lpd"), n.iter=1000000, thin=1000)
}
)

#assign posterior summary statistics to an object called stats
save(codaCP, file = "output/jm_cp_lpd.Rdata")
#statsCP <-summary(codaCP2)$stat

########################################
### NH model

#initials
initsNH=list(list(log.alpha = runif(n=Ndog, min = -1.5, max = -0.5), log.beta = runif(n=Ndog,min = -1.5, max = -0.5)),
               list(log.alpha = runif(n=Ndog, min = -2.5, max = -1.5), log.beta = runif(n=Ndog,min = -2.5, max = -1.5)),
               list(log.alpha = runif(n=Ndog, min = -4, max = -3), log.beta = runif(n=Ndog,min = -1, max = -0.1)))

# initialize model
jmNH <-jags.model("R_model_dogs_NH.R", data=dat1, inits=initsNH, n.chains = 3)
# initial coda samples 
#update the model for a sufficient number of iterations to reach convergence


#rerun coda.samples (will need to monitor lpd)
runtime.nh <- system.time(
{
#update(jmNH, 8000)
update(jmNH, 10000)
codaNH<-coda.samples(jmNH, variable.names = c("deviance", "alpha", "beta", "alpha.diff","beta.diff",
                                               "mu.alpha","mu.beta","tau.log.alpha","tau.log.beta", "lpd"), n.iter=1000000, thin = 1000)
}
save(codaNH, file = "output/jm_nh_lpd.Rdata")


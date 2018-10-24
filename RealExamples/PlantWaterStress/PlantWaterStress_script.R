###############################################################################
# Load the necessary libraries
# REMEMBER TO SET THE WORKING DIRECTORY, the folder containing 
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
if("mcmcplots" %in% installed.packages()[,1] == FALSE){
	install.packages("mcmcplots", repos="https://cran.cnr.berkeley.edu/")
	library(mcmcplots)
} else{
	library(mcmcplots)
}

###############################################################################
# Variables controlling model runs
###############################################################################

# Set variables controlling the number of update iterations to update the jags 
# model when initializing it, the length of the burnin for the MCMC chains, and how
# many iterations to run (update) the jags model for.
# n.chains - The number of chains in the MCMC
# n.burnin - The length of burnin to be used
# adaptits - The number of iterations for adaptation (tuning) the MCMC
# hb\cp\nhits - The number of iterations the hb, cp, and nh models will be run
# hb\cp\nhthin - The thinning interval to be applied when running the hb, cp, and nh models.
n.chains <- 3
n.burnin <- 100
adaptits <- 1000

hbn.burnin <- 100
hbits <- 200000
hbthin <- 100

cpn.burnin <- 100
cpits <- 200000
cpthin <- 100

nhn.burnin <- 100
nhits <- 200000
nhthin <- 100


###############################################################################
# Load data and run the HB, CP, and NH models
###############################################################################


# Load data and reshape it for running the model
load("data/wpData.Rdata")

attach(baseline)
#str(baseline)
#lowest value
#min(baseline$wp)

# Set means for covariate centering of covariates in regression model.
meanmaxD<-5
meanswc30<-10
meanswc60<-13

# Number of observations
N<-nrow(baseline)
# Inputs (data) reqired for setting Dirichlet priors for importance weights in models.
alphaA <- c(1,1,1,1,1,1,1)
alphaB <- c(1,1,1,1,1,1,1)
alphaC <- c(1,1,1,1,1,1,1)

# Create list of data objects for models;
datalist<-list(wp=wp, shrubID=shrubID2,
               maxD=as.matrix(baseline[,match("maxD1",colnames(baseline)):match("maxD7", colnames(baseline))]), 
               swc30=as.matrix(baseline[,match("swc30_1",colnames(baseline)):match("swc30_7", colnames(baseline))]), 
               swc60=as.matrix(baseline[,match("swc60_1",colnames(baseline)):match("swc60_7", colnames(baseline))]),
               meanmaxD=meanmaxD, meanswc30=meanswc30, 
               meanswc60=meanswc60,
               N=N, alphaA=alphaA, alphaB=alphaB, alphaC=alphaC,
               Ab=rep(10,7))


###############################################################################
# Hierarchical Bayesian (HB) Model
###############################################################################
# Load initial values 
load("inits/inits_baseline_hb.Rdata")
load.module('dic')

### For all models, need to monitor log-pointwise density (lpd), to compute PPI. 
### Run HB model

# Initialized jags model (adapting phase)
jm_hb <- jags.model("models/R model baseline HB.R", data=datalist, 
                    inits=saved.state.hb[[2]], n.chains=n.chains, n.adapt = adaptits)
# Update jags model using coda.samples, and monitor quantities (parameters) of interest
codaHB <- coda.samples(jm_hb, variable.names = c("deviance", "tau", "sig", "tau.eps.b",
                                                  "wA", "wB", "wC", "lpd", "b.control", "b.watered", "b.diff",
                                                  "b", "mu.b", "deltaA", "deltaB", "deltaC", "p.b"),
                        n.iter=hbits, thin=hbthin, n.burnin=hbn.burnin)

# Save model (MCMC) output
save(codaHB, file = "output/coda_hb_baseline.Rdata")
# Visualize MCMC output:
mcmcplot(codaHB, parms=c("deviance", "tau", "sig", "tau.eps.b",
                         "wA", "wB", "wC",
                         "b", "mu.b"))

# Compute convergence diagnostics:
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(codaHB[[1]]))
sel <- c(sel, grep(pattern="log_lik", x=colnames(codaHB[[1]])))
for(i in 1:3){
    gelcoda[[i]] <- as.mcmc(codaHB[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

###############################################################################
# Complete Pooling (CP) Model
###############################################################################

# Load initial values 
load("inits/inits_baseline_cp.Rdata")

# Remove the Ab variable from the datalist as it is only needed for the HB model.   
if(length(grep(names(datalist), pattern="Ab")) > 0){           
    datalist <- datalist[-grep(names(datalist), pattern="Ab")]               
}

# For all models, need to monitor log-pointwise density (lpd), to compute PPI. 
# Run CP model   

# Initialized jags model (adapting phase)
jm_cp <- jags.model("models/R model baseline CP.R", data=datalist, 
                    inits=saved.state.cp[[2]], n.chains=n.chains, n.adapt = adaptits)

# Update jags model using coda.samples, and monitor quantities (parameters) of interest
# Parameters are monitored to evalute convergence, but we only need lpd for computing PPI.
codaCP <- coda.samples(jm_cp, variable.names = c("deviance", "tau", "sig",
                                                  "wA", "wB", "wC", "lpd",
                                                  "b", "deltaA", "deltaB", "deltaC"),
                        n.iter=cpits, thin=cpthin, cpn.burnin)

# Save model (MCMC) output
save(codaCP, file = "output/coda_cp_baseline.Rdata")
# Visualize MCMC output:
mcmcplot(codaCP, parms=c("deviance", "tau", "sig", "tau.eps.b",
                         "wA", "wB", "wC",
                         "b", "mu.b"))
# Compute convergence diagnostics:
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(codaCP[[1]]))
sel <- c(sel, grep(pattern="log_lik", x=colnames(codaCP[[1]])))
for(i in 1:3){
    gelcoda[[i]] <- as.mcmc(codaCP[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

###############################################################################
# Non-Hierarchical (nh) Model
###############################################################################
#load("inits_baseline.r")
load("inits/inits_baseline_nh.Rdata")

# Initialized jags model (adapting phase)
jm_nh <- jags.model("models/R model baseline NH.R", data=datalist, 
                    inits=saved.state.nh[[2]], n.chains=3, n.adapt = adaptits)

# Update jags model using coda.samples, and monitor quantities (parameters) of interest.
# Parameters are monitored to evalute convergence, but we only need lpd for computing PPI.
codaNH <- coda.samples(jm_nh, variable.names = c("deviance", "tau", "sig",
                                                  "wA", "wB", "wC", "lpd",
                                                  "b", "deltaA", "deltaB", "deltaC"),
                        n.iter=nhits, thin=nhthin, nhn.burnin)

# Save model (MCMC) output
save(codaNH, file = "output/coda_nh_baseline.Rdata")
# Visualize MCMC output:
mcmcplot(codaNH, parms=c("deviance", "tau", "sig", "tau.eps.b",
                         "wA", "wB", "wC",
                         "b", "mu.b"))
# Compute convergence diagnostics:
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(codaNH[[1]]))
for(i in 1:3){
    gelcoda[[i]] <- as.mcmc(codaNH[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

# Check convergence
gel.nh <- gelman.diag(codaNH, multivariate = F)
gel.nh$psrf[match("deviance", row.names(gel.nh$psrf)),]


###############################################################################
# Compute PPI
###############################################################################

# Compute PPI given the lpd values (coda) from the HB, CP, and NH models:
source("source/PPI.R")
ppi.data<-PPI(codahb = codaHB, codanh = codaNH, codacp = codaCP, name="lpd")

# Uncomment the below line to save the workspace
# save.image(file='PlantWaterStress_workspace.RData')
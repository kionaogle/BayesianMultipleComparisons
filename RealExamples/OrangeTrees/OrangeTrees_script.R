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

# Set variables controlling the number of update iterations to update the Jags 
# model when initializing it, the length of the burnin for the MCMC chains, and how
# many iterations to run (update) the Jags model for.
# n.chains - The number of chains in the MCMC
# n.burnin - The length of burnin to be used
# adaptits - The number of iterations for adaptation (tuning) the MCMC
# hb\cp\nhits - The number of iterations the model will be ran
# hb\cp\nhthin - The thinning interval to be applied when running the model.
n.chains <- 3
n.burnin <- 100
adaptits <- 1000

hbn.burnin <- 10000
hbits <- 2000000
hbthin <- 1000

cpn.burnin <- 10000
cpits <- 2000000
cpthin <- 1000

nhn.burnin <- 10000
nhits <- 2000000
nhthin <- 1000


###############################################################################
# Load data and run the HB, CP, and NH models
###############################################################################

indat <- read.csv("data/treegrowth.csv")
attach(indat)

dat <- list(Y=Y, Tree=Tree, Time=Time, N=nrow(indat), Ntree=length(unique(Tree)), Nparms=3)

# Number of trees and number of pairwise comparisons:
Ntree = dat$Ntree
Ncomp = Ntree*(Ntree-1)/2

# Indices for tree-tree pairwise comparisons
pairID = matrix(data =NA, nrow = Ncomp, ncol=2)
ii = 0
for(s1 in 1:(Ntree-1)){
  for(s2 in (s1+1):Ntree){
    ii = ii +1
    pairID[ii,1] = s1
    pairID[ii,2] = s2
  }
}

# Create data list for supplying to Jags:
# For HB model:
data = list(Y=Y, Tree=Tree, Time=Time, N=nrow(indat), Ntree=length(unique(Tree)), Nparms=3, Ncomp=Ncomp, pairID = pairID)
# For CP and NH models:
data.cp = list(Y=Y, Tree=Tree, Time=Time, N=nrow(indat), Ntree=length(unique(Tree)), Nparms=3)


###############################################################################
# Initialize the Jags models for all three model variants (HB, CP, and NH)

# Hierarchical model:
# Set initials:
inits.hb = list(list(log.mu.b=log(c(190,40, 20)),tau = 0.005, tau.log.b = 1/(c(0.1,0.5,2)^2)),
                list(log.mu.b=log(c(150,30, 60)),tau = 0.05, tau.log.b = 1/(c(0.5,0.2,0.5)^2)),
                list(log.mu.b=log(c(195,35, 40)),tau = 0.015, tau.log.b = 1/(c(1,1,1.5))^2))
# Initialize model:
jm.hb <- jags.model(file = "models/R model trees HB.R", data = data, n.chains = n.chains, n.adapt = adaptits, inits=inits.hb)

# Complete pooling model:
# Set initials:
inits.cp = list(list(log.mu.b=log(c(190,40, 20)),tau = 0.005),
                list(log.mu.b=log(c(150,30, 60)),tau = 0.05),
                list(log.mu.b=log(c(195,35, 40)),tau = 0.015))
# Initialize model:
jm.cp <- jags.model(file = "models/R model trees CP.R", data = data.cp, n.chains = n.chains, n.adapt = adaptits, inits=inits.cp)

# Non-hierarchical model:
# Load initial values (will load "inits")
load("inits/nh_inits.Rdata")                
# Initialize model:
jm.nh <- jags.model(file = "models/R model trees NH.R", data = data.cp, n.chains = n.chains, n.adapt = adaptits, inits=inits)


################################## Set Monitors ###############################
# Get variable names for monitoring quantities (parameters, lpd, etc.)
jm.hb.names <- c(variable.names(jm.hb), "deviance")
jm.cp.names <- c(variable.names(jm.cp), "deviance")
jm.nh.names <- c(variable.names(jm.nh), "deviance")

load.module("dic")


##### Run / update Hierarchical model #####
update(jm.hb, n.iter=hbn.burnin)
codaHB =coda.samples(jm.hb,variable.names=jm.hb.names,
                     n.iter=hbits, thin=hbthin)

save(codaHB, file = "output/jm_hb.Rdata")


##### Run / update Complete pooling model #####
update(jm.cp, n.iter=cpn.burnin)
codaCP =coda.samples(jm.cp,variable.names=jm.cp.names,
                     n.iter=cpits, thin=cpthin)

save(codaCP, file = "output/jm_cp.Rdata")	


##### Run / update Non-hierarchical model #####
update(jm.nh, n.iter=nhn.burnin)
codaNH =coda.samples(jm.nh,variable.names=jm.nh.names,
                     n.iter=nhits, thin=nhthin)

save(codaNH, file = "output/jm_nh.Rdata")			  


# Plot some output to visulaize MCMC results and to 
# qualitatively check for convergence (see Plant_water_stress 
# example for additional code for evaluating convergence):
mcmcplot(codaHB, parms=c("deviance", "b", "mu.b", "b.diff","sig","sig.log.b"))
mcmcplot(codaCP, parms=c("deviance", "mu.b", "sig"))
mcmcplot(codaNH, parms=c("deviance", "b", "sig"))


###############################################################################
# Compute PPI
###############################################################################

# Compute PPI given the lpd values (stored in coda) from the HB, CP, and NH models:
source("source/PPI.R")
ppi.data <- PPI(codahb = codaHB, codanh = codaNH, codacp = codaCP, name="lpd")

# Uncomment the below line to save the workspace
# save.image(file='OrangeTrees_workspace.RData')
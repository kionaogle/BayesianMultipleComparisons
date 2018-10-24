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
n.burnin <- 500
adaptits <- 10000

hbn.burnin <- 500
hbits <- 50000
hbthin <- 25

cpn.burnin <- 500
cpits <- 1000000
cpthin <- 500

nhn.burnin <- 500
nhits <- 1000000
nhthin <- 500


###############################################################################
# Load data and run the HB, CP, and NH models
###############################################################################

# Load data
d <- read.csv("data/dogshock_matrix.csv", header = F)
Y = as.matrix(d)
#dim(Y)

# Number of dogs and number of pairwise comparisons:
Ndog = nrow(d)
Ncomp = Ndog*(Ndog-1)/2

# Indices for dog pairwise comparisons
pairID = matrix(data =NA, nrow = Ncomp, ncol=2)
ii = 0
for(s1 in 1:(Ndog-1)){
  for(s2 in (s1+1):Ndog){
    ii = ii +1
    pairID[ii,1] = s1 # first dog in pair
    pairID[ii,2] = s2 # second dog in pair
  }
}

# Create the xa and xs matrices of covariates (total past successes & total past failures)
xa <- matrix(NA, nrow=30, ncol=25)
xs <- matrix(NA, nrow=30, ncol=25)
for(i in 1:30){
  # Set starting values for the past number of successes (avoiding shock)
  xa[i, 1] <- 0
  # For the past number of failures (getting shocked)
  xs[i, 1] <- 0
    for(t in 2:25){
    # Sum number of past successes (avoiding shock) prior to trial t
    xa[i, t] <- sum(Y[i, 1:(t-1)])
    # Sum number of past failures (shocked) prior to trial t
    xs[i, t] <- (t-1) - xa[i,t]
  }
}

# Create data list for supplying to Jags:
datain <- list(Y=as.matrix(d), Ndog=nrow(d), Ntrial=ncol(d), xa=xa, xs=xs, pairID = pairID, Ncomp = Ncomp)

# Initials for HB model:
initsHB <- list(list(mu.log.alpha = -1, mu.log.beta = -1, tau.log.beta = 10, tau.log.alpha = 10),
                 list(mu.log.alpha = 1, mu.log.beta = 0, tau.log.beta = 20, tau.log.alpha = 50),
                 list(mu.log.alpha = -.1, mu.log.beta = .8, tau.log.beta = 1, tau.log.alpha = 5))
                 
initsNH <- list(list(log.alpha = runif(n=Ndog, min = -1.5, max = -0.5), log.beta = runif(n=Ndog,min = -1.5, max = -0.5)),
               list(log.alpha = runif(n=Ndog, min = -2.5, max = -1.5), log.beta = runif(n=Ndog,min = -2.5, max = -1.5)),
               list(log.alpha = runif(n=Ndog, min = -4, max = -3), log.beta = runif(n=Ndog,min = -1, max = -0.1)))
               
# Initials for CP model:
initsCP <- list(list(log.alpha = -1, log.beta = -1),
               list(log.alpha = 1, log.beta = 0),
               list(log.alpha = -.1, log.beta = .8))

# Create data for the CP and NH models by removing quantities from datain that are not 
# required for these models:
cpdata <- datain[-c(6, 7)]
nhdata <- datain[-c(6, 7)]

###############################################################################
# Initialize the Jags models for all three model variants (HB, CP, and NH)

# Hierarchical model:
jm.hb <- jags.model(file = "models/R_model_dogs_HB.R", data = datain, inits = initsHB, n.chains = n.chains, n.adapt = adaptits)

# Complete pooling model:
jm.cp <- jags.model(file = "models/R_model_dogs_CP.R", data = cpdata, inits = initsCP, n.chains = n.chains, n.adapt = adaptits)

# Non-hierarchical model:
jm.nh <- jags.model(file = "models/R_model_dogs_NH.R", data = nhdata, inits = initsNH, n.chains = n.chains, n.adapt = adaptits)

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
save(codaHB, file = "output/jm_hb_.Rdata")
			  
##### Run / update Complete pooling model #####
update(jm.cp, n.iter=cpn.burnin)
codaCP =coda.samples(jm.cp,variable.names=jm.cp.names,
                  n.iter=cpits, thin=cpthin)
save(codaCP, file = "output/jm_cp_.Rdata")	

##### Run / update Non-hierarchical model #####
update(jm.nh, n.iter=nhn.burnin)
codaNH =coda.samples(jm.nh,variable.names=jm.nh.names,
                  n.iter=nhits, thin=nhthin)
save(codaNH, file = "output/jm_nh_.Rdata")

# Plot some output to visulaize MCMC results and to 
# qualitatively check for convergence (see Plant_water_stress 
# example for additional code for evaluating convergence):
mcmcplot(codaHB, parms=c("deviance", "alpha", "beta"))
mcmcplot(codaCP, parms=c("deviance", "alpha", "beta"))
mcmcplot(codaNH, parms=c("deviance", "alpha", "beta"))


###############################################################################
# Compute PPI
###############################################################################

# Compute PPI given the lpd values (stored in coda) from the HB, CP, and NH models:
source("source/PPI.R")
ppi.data <- PPI(codahb = codaHB, codanh = codaNH, codacp = codaCP, name="lpd")

# Uncomment the below line to save the workspace
# save.image(file='DogLearning_workspace.RData')

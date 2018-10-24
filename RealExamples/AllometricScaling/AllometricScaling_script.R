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

# Data collected by Chuck Price (Price et al. [2007] PNAS 104:13204-13209; 
# Price et al. [2009] Ecology Letters 12:641-651)
chuck = read.csv("data/ChuckDataMatrix.csv")

# Number of species and number of pairwise comparisons:
Nsp = 49
Ncomp = Nsp*(Nsp-1)/2

# Indices for species pairwise comparisons
pairID = matrix(data =NA, nrow = Ncomp, ncol=2)
ii = 0
for(s1 in 1:(Nsp-1)){
  for(s2 in (s1+1):Nsp){
    ii = ii +1
    pairID[ii,1] = s1 # first species in pair
    pairID[ii,2] = s2 # second species in pair
  }
}

# Create data list for supplying to Jags, for HB model:
data = list(N=1162, Nsp=Nsp, R=structure(.Data=c(1,0,0,1),.Dim=c(2,2)), pairID = pairID, Ncomp = Ncomp,
            SP = chuck$SP, LogMass = chuck$LogMass, LogLength = chuck$LogLength, LogDiameter = chuck$LogDiameter,
            nVars=2, pi=pi)

# Create data list for supplying to Jags, for CP and NH models:
data_cp = list(N=1162, Nsp=Nsp, R=structure(.Data=c(1,0,0,1),.Dim=c(2,2)),
            SP = chuck$SP, LogMass = chuck$LogMass, LogLength = chuck$LogLength, LogDiameter = chuck$LogDiameter,
            nVars=2, pi=pi)

# Get initials for supplying to Jags
# For HB model (will load chuck_inits):
source("inits/R inits chuck HB.R")
# For CP model (will load chuck_inits_CP):
source("inits/R inits chuck CP.R")
# For NH model (will load chuck_inits_NH):
load("inits/R inits chuck NH.Rdata")

###############################################################################
# Initialize the Jags models for all three model variants (HB, CP, and NH)

# Hierarchical model:
jm.hb <- jags.model(file = "models/R model HB.R", data = data, inits = chuck_inits, n.chains = n.chains, n.adapt = adaptits)

# Complete pooling model:
jm.cp <- jags.model(file = "models/R model CP.R", data = data_cp, inits = chuck_inits_CP, n.chains = n.chains, n.adapt = adaptits)

# Non-hierarchical model:
jm.nh <- jags.model(file = "models/R model NH.R", data = data_cp, inits = chuck_inits_NH, n.chains = n.chains, n.adapt = adaptits)

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

##### Load model objects to output folder #####				  
#load(file = "output/jm_hb.Rdata")
#load(file = "output/jm_cp.Rdata")
#load(file = "output/jm_nh.Rdata")
# uncomment the above three lines to load from files

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
# save.image(file='AllometricScaling_workspace.RData')

###############################################################################
# for R2Bugs packages
###############################################################################

# This was a function but that is not needed 
print(paste0("multcomp function call ", sim, " started."))

#load libraries
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
#library(coda, lib.loc = "/home/mkf58/R/x86_64-redhat-linux-gnu-library/3.5")
#library(rjags, lib.loc = "/home/mkf58/R/x86_64-redhat-linux-gnu-library/3.5")
#library(xtable, lib.loc = "/home/mkf58/R/x86_64-redhat-linux-gnu-library/3.5")
#library(coda)
#library(rjags)

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

dat1 <- list(Y=as.matrix(d), Ndog=nrow(d), Ntrial=ncol(d), xa=xa, xs=xs, pairID = pairID, Ncomp = Ncomp)
datain <- dat1
print(paste0(fileload, " loaded."))

# initials
initsHB <-list(list(mu.log.alpha = -1, mu.log.beta = -1, tau.log.beta = 10, tau.log.alpha = 10),
                 list(mu.log.alpha = 1, mu.log.beta = 0, tau.log.beta = 20, tau.log.alpha = 50),
                 list(mu.log.alpha = -.1, mu.log.beta = .8, tau.log.beta = 1, tau.log.alpha = 5))

n.chains <- 3
chuck_inits <- initsHB
###############################################################################
# Hierarchical model:
jm.hb <- jags.model(file = "models/R_model_dogs_HB.R", data = datain, inits = chuck_inits, n.chains = n.chains, n.adapt = 1000)

# Complete pooling model:
jm.cp <- jags.model(file = "models/R_model_dogs_CP.R", data = datain, inits = chuck_inits, n.chains = n.chains, n.adapt = 1000)

# Non-hierarchical model:
jm.nh <- jags.model(file = "models/R_model_dogs_NH.R", data = datain, inits = chuck_inits, n.chains = n.chains, n.adapt = 1000)

################################## Added by MKF ###############################
# 8/31/2018
# Burn In
jm.hb.names <- c(variable.names(jm.hb), "deviance")
jm.cp.names <- c(variable.names(jm.cp), "deviance")
jm.nh.names <- c(variable.names(jm.nh), "deviance")

load.module("dic")

##### Run Hierarchical model #####
n.iters = 10000
n.itersrun = 50000
thinint = 25

# For testing
#n.iters = 1000
#n.itersrun = 2000
#thinint = 1

print("run hb model")
runtime.hb <- system.time(
{
update(jm.hb, n.iter=n.iters)
codahb =coda.samples(jm.hb,variable.names=jm.hb.names,
                  n.iter=n.itersrun, thin=thinint)
}
)
filename <- paste0("output/jm_hb_.Rdata")
save(codahb, file = filename)
print(runtime.hb)
rm(list=c("codahb"))
gc()

n.iters = 10000
n.itersrun = 1000000
thinint = 500
print("run cp model")				  
##### Run Complete pooling model #####
runtime.cp <- system.time(
{
update(jm.cp, n.iter=n.iters)
codacp =coda.samples(jm.cp,variable.names=jm.cp.names,
                  n.iter=n.itersrun, thin=thinint)
}
)		
filename <- paste0("output/jm_cp_.Rdata")
save(codacp, file = filename)	
print(runtime.cp)
rm(list=c("codacp"))
gc()


n.iters = 10000
n.itersrun = 1000000
thinint = 500
print("run nh model")
##### Run Non-hierarchical model #####
runtime.nh <- system.time(
{
update(jm.nh, n.iter=n.iters)
codanh =coda.samples(jm.nh,variable.names=jm.nh.names,
                  n.iter=n.itersrun, thin=thinint)
}
)
filename <- paste0("output/jm_nh_.Rdata")	
save(codanh, file = filename)
print(runtime.nh)			  

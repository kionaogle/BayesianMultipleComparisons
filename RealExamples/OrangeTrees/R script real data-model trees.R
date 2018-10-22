setwd("Z:\\common\\Projects\\mult_comparison\\Github\\BayesianMultipleComparisons\\RealExamples\\OrangeTrees")

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)

indat<-read.csv("treegrowth.csv")
names(indat)
attach(indat)

dat<-list(Y=Y, Tree=Tree, Time=Time, N=nrow(indat), Ntree=length(unique(Tree)), Nparms=3)

Ntree = dat$Ntree
Ncomp = Ntree*(Ntree-1)/2

# Indices for species pairwise comparisons; number of comparisons = (Nsp)*(Nsp-1)/2
pairID = matrix(data =NA, nrow = Ncomp, ncol=2)
ii = 0
for(s1 in 1:(Ntree-1)){
  for(s2 in (s1+1):Ntree){
    ii = ii +1
    pairID[ii,1] = s1
    pairID[ii,2] = s2
  }
}

data = list(Y=Y, Tree=Tree, Time=Time, N=nrow(indat), Ntree=length(unique(Tree)), Nparms=3, Ncomp=Ncomp, pairID = pairID)

n.chains = 3

# Hierarchical model:
runtime.hb <- system.time(
{
jm.hb <- jags.model(file = "models/R model trees HB.R", data = data, n.chains = n.chains, n.adapt = 1000000)
# initial coda samples 
codaHB<-coda.samples(jm.hb, variable.names = c("deviance", "b", "mu.b", "b.diff","sig","sig.log.b", "p.b", "lpd"), n.iter=2000000, thin=1000)
}
)
print(runtime.hb)
save(jm.hb, file="output/jm_hb_model.Rdata")
#save(codaHB, file = "output/jm_hb.Rdata")
mcmcplot(codaHB, parms=c("b", "mu.b", "b.diff","sig","sig.log.b"))
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(codaHB[[1]]))
sel <- c(sel , grep(pattern="p.b", x=colnames(codaHB[[1]])))
for(i in 1:n.chains){
    gelcoda[[i]] <- as.mcmc(codaHB[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

# Complete pooling model:
runtime.cp <- system.time(
{
data.cp = list(Y=Y, Tree=Tree, Time=Time, N=nrow(indat), Ntree=length(unique(Tree)), Nparms=3)
jm.cp <- jags.model(file = "models/R model trees CP.R", data = data.cp, n.chains = n.chains, n.adapt = 1000000)
codaCP<-coda.samples(jm.cp, variable.names = c("deviance","mu.b", "sig", "lpd"), n.iter=2000000, thin = 1000)
}
)
print(runtime.cp)
save(codaCP, file = "output/jm_cp.Rdata")
mcmcplot(codaCP, parms=c("deviance","mu.b", "sig"))
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(codaCP[[1]]))
for(i in 1:n.chains){
    gelcoda[[i]] <- as.mcmc(codaCP[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)

# Non-hierarchical model:
#load("nh_inits_mkf_2018_10_05.Rdata")
load("output/jm_nh_model.Rdata")
#inits <- list(list(tau=0.035465898, log.b=log(inits[[2]][[1]]$b)), 
#                list(tau=0.026962826, log.b=log(inits[[2]][[2]]$b)), 
#                list(tau=0.030778701, log.b=log(inits[[2]][[3]]$b)))
#save(inits, file="nh_inits_fromhb_mkf_2018_10_05.Rdata")
load("nh_inits_fromhb_mkf_2018_10_05.Rdata")                
runtime.nh <- system.time(
{
jm.nh <- jags.model(file = "models/R model trees NH.R", data = data.cp, n.chains = n.chains, n.adapt = 1000, inits=inits)
codaNH<-coda.samples(jm.nh, variable.names = c("deviance","b", "sig", "lpd", "tau", "mu"),
                    n.iter=110000, thin=50, , n.burnin=10000)
}
)	
print(runtime.nh)
save(jm.nh, file="output/jm_nh3_model.Rdata")
save(codaNH, file = "output/jm_nh3.Rdata")	
mcmcplot(codaNH, parms=c("deviance","b", "sig"))
gelcoda <- list()
sel <- grep(pattern="lpd", x=colnames(codaNH[[1]]))
sel <- c(sel, grep(pattern="mu", x=colnames(codaNH[[1]])))
sel <- c(sel, grep(pattern="tau", x=colnames(codaNH[[1]])))
for(i in 1:n.chains){
    gelcoda[[i]] <- as.mcmc(codaNH[[i]][,-sel])
}
gelman.diag(x=gelcoda, multivariate=FALSE)
coda.fast(chains=3, thin=1, coda = gelcoda)
inits <- initfind(mcmcin=codaHB)
inituse <- removevars(inits, c(1,2,3,4))
inituse[[2]][[2]]$tau  <- (inits[[2]][[1]]$tau + inits[[2]][[3]]$tau) / 2
save(inituse, file="nh_inits_mkf_2018_10_05.Rdata")




# Target, hierarchical model for orange tree growth example
model{
  for(i in 1:N){
    # Normal likelihood for circumference data (Y)
    Y[i] ~ dnorm(mu[i], tau)
    # log pointwise density
    lpd[i] <- logdensity.norm(Y[i], mu[i], tau)
    # let b[1,i] = Ymax[i] (maximum circumference)
    # let b[2,i] = Ymin[i] (circumference at time of planting)
    # let b[3,i] = Ghalf[i] (growth rate at half max circumference)
    # Non-linear, logistic growth curve for predicted (mean) circumference:
    mu[i] <- b[1,Tree[i]]/(1+exp(-4*b[3,Tree[i]]/b[1,Tree[i]]*Time[i])*
                             (b[1,Tree[i]]/b[2,Tree[i]]-1))
  }
  # Hierachical priors for tree-level growth parameters
  for(p in 1: Nparms){
    for(t in 1:Ntree){
      # Log scale
      log.b[p,t] ~ dnorm(log.mu.b[p], tau.log.b[p])
      # Regular scale
      b[p,t] <- exp(log.b[p,t])
    }
  }
  # Conjugate, relatively non-informative hyperpriors for population-
  # level parameters:
  for(i in 1:Nparms){
    # Log scale
    log.mu.b[i] ~ dnorm(0, 0.0001)
    # Regular scale
    mu.b[i] <- exp(log.mu.b[i])
    # Precision and standard deviation describing variation among trees
    tau.log.b[i] ~ dgamma(0.01, 0.01)
    sig.log.b[i] <- pow(tau.log.b[i], -0.5)
  }
  
  # Conjugate, relativley non-informative gamma prior for tau (precition),
  # calculate sig (standard deviation)
  tau ~ dgamma(0.01, 0.01)
  sig <- pow(tau, -0.5)
  
  # Pairwise comparisons of tree growth parameters among individual trees:
  for(j in 1:Ncomp){
    for(p in 1:Nparms){
      b.diff[j,p] <- b[p,pairID[j,1]] - b[p,pairID[j,2]]
      # Posterior means of the quantities below give Bayesian 
      # p-values for each pairwise comparison
      p.b[j,p] <- step(b.diff[j,p])
    }
  }
}